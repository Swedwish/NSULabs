package ru.nsu.ccfit.saburov.view;

public interface View {

    void updateGameTextField();

    void help();

    void wrong();

    void loss();

    void victory();

    void greetings();

    String[] makeMove();
}
