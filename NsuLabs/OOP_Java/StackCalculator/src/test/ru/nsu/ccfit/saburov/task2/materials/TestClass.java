package test.ru.nsu.ccfit.saburov.task2.materials;

import org.junit.jupiter.api.*;
import ru.nsu.ccfit.saburov.task2.materials.*;

import static org.junit.jupiter.api.Assertions.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Stack;
import java.util.Vector;

class TestClass {

    void pop(Context context, InputStreamReader reader) {
        (new Pop()).makeOperation(context, reader);
    }

    void push(Context context, InputStreamReader reader) {
        (new Push()).makeOperation(context, reader);
    }

    void define(Context context, InputStreamReader reader) {
        (new Define()).makeOperation(context, reader);
    }

    void sqrt(Context context, InputStreamReader reader) {
        (new Sqrt()).makeOperation(context, reader);
    }

    void math(Context context,InputStreamReader reader) {
        (new Maths()).makeOperation(context, reader);
    }

    @Test
    void PopAndPushTest() throws FileNotFoundException {
        Context context = new Context(new Stack<>(),new HashMap<>(),new Vector<>());

        //context.argsAdd("3.0");
        var input = new FileInputStream("src/test/ru/nsu/ccfit/saburov/task2/materials/PopAndPushTest.txt.");
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 3.0);
        context.argsClear();

        //context.argsAdd("5");
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 5.0);
        context.argsClear();

        pop(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 3.0);

        pop(context, new InputStreamReader(input));
        assertTrue(context.dataEmpty());
    }

    @Test
    void DefineTest() throws FileNotFoundException {
        Context context = new Context(new Stack<>(),new HashMap<>(),new Vector<>());

        context.argsAdd("a");
        context.argsAdd("4");
        var input = new FileInputStream("src/test/ru/nsu/ccfit/saburov/task2/materials/DefineTest.txt");
        define(context, new InputStreamReader(input));
        context.argsRemove(1);
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 4.0);
        context.argsClear();
        context.dataClear();

        context.argsAdd("5.0");
        context.argsAdd("4");
        define(context, new InputStreamReader(input));
        context.argsRemove(1);
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 4.0);
        context.argsClear();
        context.dataClear();

        context.argsAdd("+5");
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), 5.0);
        context.argsClear();
        context.dataClear();

        context.argsAdd("-5");
        push(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(), -5.0);
        context.argsClear();
        context.dataClear();
    }

    @Test
    void MathTest() throws FileNotFoundException {
        Context context = new Context(new Stack<>(),new HashMap<>(),new Vector<>());
        var input = new FileInputStream("MathTest.txt");
        context.argsAdd("25");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"5.0");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"+");
        math(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(),30);
        context.argsClear();

        context.argsAdd("25");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"5.0");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"-");
        math(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(),20);
        context.argsClear();
        
        context.argsAdd("25");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"5.0");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"/");
        math(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(),5);
        context.argsClear();

        context.argsAdd("25");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"5.0");
        push(context, new InputStreamReader(input));
        context.argsSet(0,"*");
        math(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(),125);
        context.argsClear();

        context.argsAdd("25");
        push(context, new InputStreamReader(input));
        context.argsClear();
        sqrt(context, new InputStreamReader(input));
        assertEquals(context.dataPeek(),5);
    }
}