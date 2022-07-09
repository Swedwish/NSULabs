#ifndef LAB1_HEADER_H
#define LAB1_HEADER_H
#pragma once

#include <algorithm>
#include "TritSet.h"

typedef unsigned int uint;


using namespace std;

int bit(uint b, uint bitNumber) {
    return (b & (1 << bitNumber)) >> bitNumber;
}

TritSet operator&(const TritSet &first, const TritSet &second) {
    TritSet temp;
    temp = first;
    temp &= second;
    return temp;
}

TritSet operator|(const TritSet &first, const TritSet &second){
    TritSet temp;
    temp = first;
    temp |= second;
    return temp;
}

static int getTrit(uint b, uint bitNumber) {
    if (bit(b, bitNumber) == 1) {
        return True;
    } else if (bit(b, bitNumber + 1) == 1) {
        return False;
    } else {
        return Unknown;
    }
}


TritSet::reference::reference() {
    index = 0;
}

TritSet::reference& TritSet::reference::operator=(Trit value) {
    uint bit_index = index * 2 / (8 * sizeof(uint));
    if (index >= last->capacity()) {
        if (value == True || value == False) {
            last->resize(index + 1);
        } else
            return *this;
    }
    if (value == False) {
        last->data[bit_index] |= 1 << (index % 16 * 2 + 1);//second bit index
    } else {
        last->data[bit_index] &= ~(1 << (index % 16 * 2 + 1));//second bit index
    }
    if (value == True) {
        last->data[bit_index] |= 1 << (index % 16 * 2); //first bit index
    } else {
        last->data[bit_index] &= ~(1 << (index % 16 * 2)); //first bit index
    }
    return *this;
}

bool TritSet::reference::operator==(Trit A) {
    if (last->length > index) {
        return (getTrit(last->data[index * 2 / (8 * sizeof(uint))], index * 2 % (8 * sizeof(uint))) == A);
    } else {
        return (A == Unknown);
    }
}

uint *TritSet::show_data() {
    return this->data;
}


uint TritSet::capacity() const {
    return length;
}

void TritSet::shrink() {
    int last_one = 0, flag = 0;
    for (int i = 0; i < capacity() * 2; i += 2) {
        uint byte_index = ((i) / (8 * sizeof(uint)));
        uint bit_index = ((i) % (8 * sizeof(uint)));
        if (bit(data[byte_index], bit_index) == 1 || bit(data[byte_index], bit_index + 1) == 1) {
            last_one = i / 2;
            flag = 1;
        }
    }
    if (flag == 1) {
        resize(last_one + 1);
    } else {
        resize(0);
    }
}

void TritSet::resize(int n_len) {
    TritSet Trit_ban(n_len);
    Trit_ban.data = new uint[n_len * 2 / (8 * sizeof(uint))];
    int temp = length != 0 ? 1 : 0;
    int kekwwww = n_len * 2 % (8 * sizeof(uint)) == 0 ? 0 : 1;
    temp += length * 2 / (8 * sizeof(uint));
    if (temp > n_len * 2 / (8 * sizeof(uint)) + kekwwww) {
        temp = n_len * 2 / (8 * sizeof(uint)) + kekwwww;
    }
    //temp = min(temp, n_len * 2 / (8 * sizeof(uint))+kekwwww);
    for (uint i = 0; i < (temp); i++) {
        Trit_ban.data[i] = data[i];
    }
    for (uint i = length * 2; i < n_len * 2; i++) {
        Trit_ban.data[i / (sizeof(uint) * 8)] &= ~(1 << i % (sizeof(uint) * 8));
    }
    //delete [] data;
    data = Trit_ban.data;
    length = Trit_ban.length;
}

TritSet::TritSet(uint len) {
    length = len;
    if (len == 0){
        data = new uint[1];
        data[0] = 0;
        return;
    }
    data = new uint[len * 2 / (8 * sizeof(uint))];
    //data = (uint*) malloc(len*2/(8*sizeof (uint)));
    int temp = len*2%8*sizeof(uint) == 0? 0:1;
    for (uint i = 0; i < (len * 2 / (8 * sizeof(uint))+temp); i++) {
        data[i] = 0;
    }
}

TritSet::TritSet() {
    length = 0;
    data = nullptr;
}

TritSet::TritSet(const TritSet& set):data(set.data),length(set.length){}

TritSet::reference kekw;

TritSet::reference &TritSet::operator[](uint index) {
    kekw.last = this;
    kekw.index = index;
    return kekw;
}

TritSet& TritSet::operator=(const TritSet& tritBan) {
    if (tritBan.length == 0){
        data = new uint[1];
        data[0] = 0;
        return *this;
    }
    uint* old_data = data;
    data = new uint[tritBan.length * 2 / (8 * sizeof(uint))];
    int temp = tritBan.length * 2 % 8 * sizeof(uint) == 0 ? 0 : 1;
    for (uint i = 0; i < (length * 2 / (8 * sizeof(uint)) + temp); i++) {
        data[i] = 0;
    }
    temp = tritBan.length * 2 % 8 * sizeof(uint) == 0 ? 0 : 1;
    for (uint i = 0; i < (tritBan.length * 2 / (8 * sizeof(uint)) + temp); i++) {
        this->data[i] = tritBan.data[i];
    }
    this->length = tritBan.length;
    delete[] old_data;
    return *this;
}

int TritSet::logical_length() {
    int last_one = 0;
    for (int i = 0; i < capacity() * 2; i += 2) {
        uint byte_index = ((i) / (8 * sizeof(uint)));
        uint bit_index = ((i) % (8 * sizeof(uint)));
        if (bit(data[byte_index], bit_index) == 1 || bit(data[byte_index], bit_index + 1) == 1) {
            last_one = i / 2;
        }
    }
    return last_one + 1;
}

void TritSet::trim(int lastIndex) {
    resize(lastIndex);
}

int TritSet::cardinality(const Trit value) {
    shrink();
    int c = 0;
    for (int i = 0; i < capacity(); i++) {
        uint byte_index = ((i * 2) / (8 * sizeof(uint)));
        uint bit_index = ((i * 2) % (8 * sizeof(uint)));
        if (getTrit(data[byte_index], bit_index) == value) {
            c++;
        }
    }
    return c;
}

void TritSet::operator&=(const TritSet &second) {
    uint maxlength;
    if (this->length > second.length) {
        maxlength = this->length;
    } else {
        maxlength = second.length;
    }
    TritSet temp(maxlength);

    for (uint i = 0; i < (maxlength); i++) {
        uint byte_index = ((i * 2) / (8 * sizeof(uint)));
        uint bit_index = ((i * 2) % (8 * sizeof(uint)));
        int ftrit, strit;
        if (this->capacity() <= i) ftrit = Unknown;
        else ftrit = getTrit(this->data[byte_index], bit_index);
        if (second.capacity() <= i) strit = Unknown;
        else strit = getTrit(second.data[byte_index], bit_index);

        if ((ftrit == False) || (strit == False)) {
            temp[i] = False;
        } else if (ftrit == Unknown || strit == Unknown) {
            temp[i] = Unknown;
        } else {
            temp[i] = True;
        }
    }
    *this = temp;
}

void TritSet::operator|=(const TritSet &second) {
    uint maxlength;
    if (this->length > second.length) {
        maxlength = this->length;
    } else {
        maxlength = second.length;
    }
    TritSet temp(maxlength);
    for (uint i = 0; i < (maxlength); i++) {
        uint byte_index = ((i * 2) / (8 * sizeof(uint)));
        uint bit_index = ((i * 2) % (8 * sizeof(uint)));

        int ftrit, strit;
        if (this->capacity() <= i)
            ftrit = Unknown;
        else
            ftrit = getTrit(this->data[byte_index], bit_index);
        if (second.capacity() <= i)
            strit = Unknown;
        else
            strit = getTrit(second.data[byte_index], bit_index);
        if ((ftrit == True) || (strit == True)) {
            temp[i] = True;
        } else if (ftrit == Unknown || strit == Unknown) {
            temp[i] = Unknown;
        } else {
            temp[i] = False;
        }
    }
    *this = temp;
}

TritSet TritSet::operator!() {
    TritSet temp(this->length);
    for (uint i = 0; i < (this->length); i++) {
        uint byte_index = ((i * 2) / (8 * sizeof(uint)));
        uint bit_index = ((i * 2) % (8 * sizeof(uint)));
        int firstTrit;
        if (this->capacity() < i) firstTrit = Unknown;
        else firstTrit = getTrit(this->data[byte_index], bit_index);
        if (firstTrit == False) {
            temp[i] = True;
        } else if (firstTrit == Unknown) {
            temp[i] = Unknown;
        } else {
            temp[i] = False;
        }
    }
    return temp;
}

#endif //LAB1_HEADER_H
