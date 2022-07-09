#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch2.hpp" //Хэдер с сайта
#include "TritSet.cpp" //Тут лежит тритсет и всякие штуки
TEST_CASE("Default_tests", "[TritSets]") {
    //резерв памяти для хранения 1000 тритов
    TritSet set(1000);
    // length of internal array
    int allocLength = set.capacity();
    REQUIRE(allocLength >= (1000 * 2 / 8 / sizeof(uint)));
    // 1000*2 - min bits count
    // 1000*2 / 8 - min bytes count
    // 1000*2 / 8 / sizeof(uint) - min uint[] size

    //не выделяет никакой памяти
    set[1000000000] = Unknown;
    REQUIRE(allocLength == set.capacity());

    // false, but no exception or memory allocation
    if (set[2000000000] == True) {}
    REQUIRE(allocLength == set.capacity());

    //выделение памяти
    set[1000000000] = True;
    REQUIRE(allocLength < set.capacity());

    //no memory operations
    allocLength = set.capacity();
    set[1000000000] = Unknown;
    set[1000000] = False;
    REQUIRE(allocLength == set.capacity());

    //освобождение памяти до начального значения или
    //до значения необходимого для хранения последнего установленного трита
    //в данном случае для трита 1000’000
    set.shrink();
    REQUIRE(allocLength > set.capacity());

}

TEST_CASE("Operators"){    //data - поле данных структуры. В моем случае это массив uint
    TritSet A(1000);
    TritSet B(2000);
    TritSet C = A & B;
    REQUIRE(C.capacity() == B.capacity());
    A[1] = True;
    B[1] = False;
    REQUIRE((A & B).show_data()[0] == B.show_data()[0]);
    REQUIRE((A | B).show_data()[0] == A.show_data()[0]);
    REQUIRE((!A).show_data()[0] == B.show_data()[0]);
    A[0] = True; A[2] = True; A[3] = False;A[4] = False; A[5] = False;
    B[0] = Unknown; B[2] = True; B[3] = Unknown; B[4] = False; B[5] = True; B[6] = Unknown; B[7] = False; B[8] = True;
    C = A|B;
    REQUIRE((C).show_data()[0] == 67093);  //==TUUTFUTTT
    C = A&B;
    REQUIRE((C).show_data()[0] == 35480);  //==UFUFFFTFU
    TritSet D;
    TritSet E(3);
    E[0] = True; E[2] = False;
    TritSet G = E|D;
    REQUIRE(G.show_data()[0] == 1);  //==UUT
    G = E&D;
    REQUIRE(G.show_data()[0] == 32); //==FUU
}

TEST_CASE("SHRINK"){
    TritSet A(10);
    A.shrink();
    REQUIRE(A.capacity() == 0);
    A[0] = True;
    REQUIRE(A.capacity() == 1);
    A[50] = True;
    A.resize(200);
    A.shrink();
    REQUIRE(A.capacity()==51);
    A.shrink();
    REQUIRE(A.capacity()==51);
}


TEST_CASE("CARDINALITY"){
    TritSet A;
    A[5] = True;
    A[7] = False;
    int C = A.cardinality(True);
    REQUIRE(C == 1);
    C = A.cardinality(False);
    REQUIRE(C == 1);
    C = A.cardinality(Unknown);
    REQUIRE(C == 6);
    A[5] = False;
    A[9] = False;
    C = A.cardinality(True);
    REQUIRE(C == 0);
    C = A.cardinality(False);
    REQUIRE(C == 3);
    C = A.cardinality(Unknown);
    REQUIRE(C == 7);
}

TEST_CASE("TRIM"){
    TritSet B(100);
    B[10] = True;
    B[5] = False;
    B[4] = True;
    B.trim(5);
    bool A = (B[4] == True);
    REQUIRE(A);
    A = B[10] == Unknown;
    REQUIRE(A);
    A = B[5] == Unknown;
    REQUIRE(A);
    A = B[4] == True;
    REQUIRE(A);
    REQUIRE(B.capacity() == 5); //capacity в тритах, насколько я понял
}

TEST_CASE("LOGICAL LENGTH"){
    TritSet A(5);
    //print(const TritSet& A);
    A[4] = True;
    A[3] = False;
    A[6] = Unknown;
    REQUIRE(A.logical_length() == 5);
}
