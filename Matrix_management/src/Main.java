import java.io.File;

public class Main {


    /**
     * Starts the program
     * @param args list of tissue files with CpG frequency
     */
    public static void main(String[] args) {
        if (argsCheck(args)) {
            matrix_average MA = new matrix_average(args);
            MA.getAverage();
        }
    }




    /**
     * Arguments file existence check
     *
     */
    private static boolean argsCheck(String[] args) {
        try {
            for (String arg : args) {
                File file = new File(arg);
                if (!file.exists()) {
                    throw new IllegalArgumentException();
                }
            }
        } catch (IllegalArgumentException existEx) {
            System.out.println("One or more of the files do not exist!");
            return false;
        }
        return true;
    }
}
