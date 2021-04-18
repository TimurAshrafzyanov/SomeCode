#include <iostream>
#include <vector>

//Реализовать хэш-таблицу
//Для разрешения коллизий использовать квадратичное пробирование

struct Point{
    std::string value;
    bool isDeleted;
    bool isEmpty;

    Point() : value(), isDeleted(false), isEmpty(true){}
};

class HashTable{
private:
    std::vector<Point> table;
    int size;
    int capacity;

    int Hash(std::string str);
    int Hash(std::string str, int i);
    void Rehash();

public:
    bool Add(std::string str);
    bool Search(std::string str);
    bool Remove(std::string str);

    HashTable(){
        size = 0;
        capacity = 8;
        table.resize(8);
    }
};

void HashTable::Rehash(){
    std::vector<Point> oldtable(capacity);
    oldtable = table;
    std::vector<Point> newtable(2 * capacity);
    table = newtable;
    capacity = capacity * 2;
    for (int i = 0; i < oldtable.size(); i++){
        if (!oldtable[i].isEmpty){
            Add(oldtable[i].value);
        }
    }
}

int HashTable::Hash(std::string str){
    int sum = 0;
    int k = 1;
    for (int i = 0; i < str.size(); i++){
        sum = (sum + k * int(str[i])) % capacity;
        k = (k * 31) % capacity;
    }
    return sum;
}

int HashTable::Hash(std::string str, int i){
    return (Hash(str) + i * (i + 1) / 2) % capacity;
}

bool HashTable::Add(std::string str){
    int i = 0;
    int ind = Hash(str, i);
    while (!table[ind].isEmpty){
        i += 1;
        if ((!table[ind].isDeleted) && (table[ind].value == str)) return false;
        ind = Hash(str, i);
    }
    table[ind].value = str;
    table[ind].isEmpty = false;
    if (table[ind].isDeleted) table[ind].isDeleted = false;
    size += 1;
    if (4 * size >= 3 * capacity) Rehash();
    return true;
}

bool HashTable::Search(std::string str){
    int i = 0;
    int ind = Hash(str, i);
    while ((!table[ind].isEmpty) || (table[ind].isDeleted)){
        if (table[ind].value == str)
            return !table[ind].isDeleted;
        i += 1;
        ind = Hash(str, i);
    }
    return false;
}

bool HashTable::Remove(std::string str){
    int i = 0;
    int ind = Hash(str, i);
    while ((!table[ind].isEmpty) || (table[ind].isDeleted)){
        if (table[ind].value == str) {
            if (table[ind].isDeleted) {
                return false;
            } else {
                table[ind].isDeleted = true;
                table[ind].isEmpty = true;
                return true;
            }
        }
        i += 1;
        ind = Hash(str, i);
    }
    return false;
}

int main() {
    HashTable hashTable;
    std::string str;
    std::vector<bool> result;
    char operation;
    bool flag;
    while (std::cin >> operation){
        std::cin >> str;
        if (operation == '+') {
            flag = hashTable.Add(str);
        } else if (operation == '?'){
            flag = hashTable.Search(str);
        } else flag = hashTable.Remove(str);
        result.push_back(flag);
    }
    for (int i = 0; i < result.size(); i++){
        if (result[i] == true) {
            std::cout << "OK" << std::endl;
        } else std::cout << "FAIL" << std::endl;
    }
    return 0;
}