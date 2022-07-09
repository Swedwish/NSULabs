package ru.nsu.ccfit.saburov;
//import org.lwjgl.*;
//import org.lwjgl.glfw.*;
//import org.lwjgl.opengl.*;
//import org.lwjgl.system.*;
//
//import java.nio.*;
//
//import static org.lwjgl.glfw.Callbacks.*;
//import static org.lwjgl.glfw.GLFW.*;
//import static org.lwjgl.opengl.GL11.*;
//import static org.lwjgl.system.MemoryStack.*;
//import static org.lwjgl.system.MemoryUtil.*;

import java.io.IOException;
import java.io.InputStreamReader;

public class Main {
    public static void main(String[] args) throws IOException {
        InputStreamReader reader = new InputStreamReader(System.in);
        System.out.println("Type \"g\" for GuiView or \"t\" for text interface.");
        int a = reader.read();
        new Controller().play(a);
    }
}