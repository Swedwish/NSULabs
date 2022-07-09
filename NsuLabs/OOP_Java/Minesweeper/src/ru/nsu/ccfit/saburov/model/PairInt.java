package ru.nsu.ccfit.saburov.model;

public class PairInt {
    int first;
    int second;

    public PairInt(int first, int second) {
        this.first = first;
        this.second = second;
    }

    Boolean isEqual(PairInt Pair2) {
        return this.first == Pair2.first && this.second == Pair2.second;
    }

    void setFirst(int value) {
        first = value;
    }

    void setSecond(int value) {
        second = value;
    }

    public int getFirst() {
        return first;
    }

    public int getSecond() {
        return second;
    }
}
