#ifndef LAB1_0_TRITSET_H
#define LAB1_0_TRITSET_H
#pragma once

#include <vector>
typedef unsigned int uint;

enum Trit {
    Unknown,
    True,
    False
};

class TritSet{
private:
    uint length;
    uint *data;

public:

    class reference{
        friend class TritSet;
    private:
        uint index;
        TritSet *last;

    public:
        reference();
        reference& operator=(Trit value);
        bool operator==(Trit A);
    };
    TritSet();
    explicit TritSet(uint len);
    TritSet(const TritSet& set);
    uint* show_data();
    uint capacity() const;
    void shrink();
    void resize(int n_len);
    reference &operator[](uint index);
    int logical_length();
    void trim(int lastIndex);
    int cardinality(Trit value);
    TritSet& operator=(const TritSet& tritBan);
    void operator&=(const TritSet& second);
    void operator|=(const TritSet& second);
    TritSet operator!();
};

TritSet operator&(const TritSet &first, const TritSet &second);
TritSet operator|(const TritSet &first, const TritSet &second);

#endif //LAB1_0_TRITSET_H
