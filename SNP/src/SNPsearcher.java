import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

// -r 500 -d 0.3 -s /cs/cbio/roei/snps151_23_10_2019.tsv -b /cs/cbio/roei/ilia/3_3_2020.tsv
public class SNPsearcher {

    // Argument flags:
    private final static String RANGE = "-r";
    private final static String DELTA = "-d";
    private final static String SNPDICT = "-s";
    private final static String BLOCK = "-b";

    // Arguments:
    private static int range;
    private static String blocksPath;
    private static String snpDictPath;
    private static float deltaFreq; // Frequency range from entropy maximum value
    private static String currentChr;

    private static ArrayList<String[]> snpBuffer_1 = new ArrayList<>();
    private static ArrayList<String[]> snpBuffer_2 = new ArrayList<>();

    private static BufferedReader snpDict;
    private static BufferedReader blocks;


//    private static SortedSet<String[]> distances = new TreeSet<>(new GenomePlaceComparator());
//    private static SortedSet<String[]> SNPs = new TreeSet<>(new GenomePlaceComparator());
    private static SortedMap<String[], String[]> SNPs = new TreeMap<>(new GenomePlaceComparator());
    private static String[] line;

    private static Set<String> visited = new HashSet<>();

    // Frequency that gives maximal entropy
    private final static float maxEntropy = 0.5f;

    // Counters:
    private static int totalCounter;
    private static int counterblocks = 0;


    private static float bestEntropy = 1;
    private static String[] bestSnp = null;

    private static boolean searchAll = false;

    /**
     * Main function that processes input arguments
     * @param args command line arguments
     */
    public static void main(String[] args) {
        if (argNumCheck(args)) {
            argSort(args);
            argsCheck();
            totalCounter = 0;
            try {
                snpDict = new BufferedReader(new FileReader(snpDictPath));
                // start and stop in block not including bounds
                blocks = new BufferedReader(new FileReader(blocksPath));
                String[] block = getNextBlock();
                block = getNextBlock();
                line = getNextLine();
                while (block != null) {
                    searcher(block);
                    block = getNextBlock();
                }
                for (String[] snp : SNPs.values()) {
                    System.out.println(snp[0] + " " + snp[1]);
                }
                snpDict.close();
                blocks.close();
            } catch (java.io.IOException exp) {
                System.out.println("Wrong file");
            }
//            save_CpG_SNP_location();
            SNPquality();
                // TODO: сравнить файлы с результатом (в этой папке и в ДАТА)
                // TODO: сделать красвую диаграмму с расстояними между CpG и SNP
        }
            System.out.println("Number of blocks: " + counterblocks);
            System.out.println("Total number of SNPs: " + totalCounter);
            System.out.println("Maximal potential score : " + getMaxScore());
    }


    private static void save_CpG_SNP_location() {
        try {
            Path resultPath = Paths.get(System.getProperty("user.dir"), "/cpg_snp.csv");
            PrintWriter out = new PrintWriter(resultPath.toString());
            for (String[] CpG : SNPs.keySet()) {
                out.print(CpG[3] + "," + CpG[1] + "," + SNPs.get(CpG)[3] + "," +  SNPs.get(CpG)[1] + "\n");
                System.out.print(CpG[3] + "," + CpG[1] + "," + SNPs.get(CpG)[3] + "," +  SNPs.get(CpG)[1] + "\n");
            }
            out.close();
        } catch (java.io.IOException exp) {
            System.out.println("CpG-SNP saving problem");
        }
    }

    private static void SNPquality() {
        try {
            Path resultPath = Paths.get(System.getProperty("user.dir"), "/SNPquality_" + range + "_" + deltaFreq + "_CTTC.csv");
            PrintWriter out = new PrintWriter(resultPath.toString());
            for (String[] CpG : SNPs.keySet()) {
                out.print(SNPs.get(CpG)[3] + "," +  SNPs.get(CpG)[7] + "," + SNPs.get(CpG)[5] + SNPs.get(CpG)[6] + "\n");
                System.out.print(SNPs.get(CpG)[3] + "," +  SNPs.get(CpG)[7] + "\n");
            }
            out.close();
        } catch (java.io.IOException exp) {
            System.out.println("CpG-SNP saving problem");
        }
    }


    /**
     * Searches SNP sites in required region
     * @param block CpG block
     */
    private static void searcher(String[] block) {
//        int counterPre = 0;
//        int counterPost = 0;
//        String[] minimalFreqSnp = null;
        int start = Integer.parseInt(block[1]);
        int end = Integer.parseInt(block[2]);
        bestEntropy = 1;
        while (line != null && !line[0].equals(currentChr)) {
            line = getNextLine();
        }
        while (line != null && Integer.parseInt(line[1]) <= (start - range)) {
            line = getNextLine();
        }
//        // Three next "while" loops can be replaced by one if it doesn't matter in which region SNP located (before block, inside, after block
//        while (line != null && Integer.parseInt(line[1]) < start) {
//            addSNP(line);
//            line = getNextLine();
//        }
//        while (line != null && Integer.parseInt(line[1]) < end) {
//            addSNP(line);
//            line = getNextLine();
//        }
//        while (line != null && Integer.parseInt(line[1]) <= end + range) {
//            addSNP(line);
//            line = getNextLine();
//        }
        while (line != null && Integer.parseInt(line[1]) < end + range) {
            if (!visited.contains(line[3])) {
                addSNP(line, block);
            }
            line = getNextLine();
        }
        if (bestEntropy != 1) {
            if (!searchAll) {
                //            SNPs.add(bestSnp);
                SNPs.put(block, bestSnp);
                visited.add(bestSnp[3]);
                totalCounter++;
                System.out.println(bestSnp[0] + " " + bestSnp[1] + " " + block[3]);
            } else {
                // TODO to print all snps that were added in this loop
            }
        }
//        totalCounter += counterPre + counterPost;
//        System.out.println("Counter pre: " + counterPre);
//        System.out.println("Counter post: " + counterPost);
    }


    /**
     * Checks if the SNP is better then previous ones in the block
     * @param line
     */
    private static void addSNP(String[] line, String[] block) {
        float currentFreq = Float.parseFloat(line[7]);
        float entropyDif = Math.abs(maxEntropy - currentFreq);
        if (searchAll && checkFreq(currentFreq)) {
            SNPs.put(block, line);
            bestEntropy = entropyDif;
            totalCounter ++;
        }
        else {
            if (checkFreq(currentFreq) && snpCheck(line[5], line[6]) && entropyDif < bestEntropy) {
                bestEntropy = entropyDif;
                bestSnp = line;
            }
        }
    }


    /**
     * Calculates maximal potential frequency of given SNP combination
     * (the lower the score, the more unique the combination)
     * @return maximal potential frequency of given SNP combination
     */
    private static float getMaxScore() {
        // DNA Profile Frequency Calculations:
        // http://www.biology.arizona.edu/human_bio/activities/blackett2/str_frequency.html
        float result = 1;
        try {
            Path resultPath = Paths.get(System.getProperty("user.dir"), "/top_SNPs.tsv");
//            Path resultPath = Paths.get("/cs/cbio/ilia/data/project_data/tissue_samples/vector_pre.csv");
            PrintWriter out = new PrintWriter(resultPath.toString());
            for (String[] line : SNPs.values()) {
                result *= Float.parseFloat(line[7]);
                for (int i = 0; i < line.length - 1; i++) {
                    out.print(line[i] + "\t");
                }
                out.print(line[line.length - 1] + "\n");
            }
            out.close();
        } catch (java.io.IOException exp) {
            System.out.println("Matrix saving problem");
        }
        return result*result; // p^2
    }


    /**
     * Checks if given frequency is in a given range from entropy maximum
     * @param freq Frequency of given SNP
     * @return True if given frequency is in a given range from entropy maximum
     *          False otherwise
     */
    private static boolean checkFreq(float freq) {
        return Math.abs(freq - maxEntropy) <= deltaFreq;
    }


    /**
     * Checks if it possible to detect SNP
     * @param org original nucleotide
     * @param alt alternative nucleotide
     * @return true if it possible to detect SNP
     */
    private static boolean snpCheck(String org, String alt) {
        String A = "A";
        String T = "T";
        String C = "C";
        String G = "G";
//        return (org.equals(A) && !alt.equals(A)) || (org.equals(T) && (alt.equals(A) || alt.equals(G))) ||
//                (org.equals(C) && (alt.equals(A)) || alt.equals(G)) || (org.equals(G) && !alt.equals(G));
        return true;
    }


    /**
     * Returns next read from the PAT file
     * @return Next read
     */
    private static String[] getNextLine() {
        String[] nextLine = null;
        if (SNPsearcher.snpBuffer_1.isEmpty()){
            try{
                String line = SNPsearcher.snpDict.readLine();
                if (line != null){
                    nextLine = line.split("\t");
                }
            } catch (java.io.IOException exp) {
                System.out.println("Get next line ERROR");
                return null;
            }
        } else {
            nextLine = SNPsearcher.snpBuffer_1.get(0);
            SNPsearcher.snpBuffer_1.remove(0);
        }
        return nextLine;
    }


    /**
     * Returns next block from the blocks file
     * @return Next block from the blocks file
     */
    private static String[] getNextBlock(){
        // TODO: Assume, that in block there are only three columns (chromosome, start and end) and they separated by \t.
        String[] line;
        String[] nextBlock = null;
        SNPsearcher.snpBuffer_1 = SNPsearcher.snpBuffer_2;
        SNPsearcher.snpBuffer_2 = new ArrayList<>();
        try{
            String block = SNPsearcher.blocks.readLine();
            if (block != null){
                line = block.split("\t");
//                currentChr = line[0].substring(3);
                currentChr = line[0];
                nextBlock = line;
                counterblocks++;
//                nextBlock = Arrays.copyOfRange(line, 1, 3);
            }
        } catch (java.io.IOException exp) {
            System.out.println("Get next block ERROR");
            return null;
        }
        return nextBlock;
    }


    /**
     * Arguments numbers legacy check
     *
     * @param args
     */
    private static boolean argNumCheck(String[] args) {
        try {
            // With result dir path
            if (args.length != 8) {
                throw new IllegalArgumentException();
//                // Without result dir path
//                if (args.length != 4) {
//                    throw new IllegalArgumentException();
//                }
            }
        } catch (IllegalArgumentException lenEx) {
            System.out.println("Illegal number of arguments!");
            return false;
        }
        return true;
    }


    /**
     * Arguments file existence check
     *
     */
    private static void argsCheck() {
        try {
            File file1 = new File(blocksPath);
            File file2 = new File(snpDictPath);
            if (!file1.exists() || !file2.exists() || deltaFreq < 0 || deltaFreq > 0.5) {
                throw new IllegalArgumentException();
            }
        } catch (IllegalArgumentException existEx) {
            System.out.println("Illegal argument / one or more of the files do not exist!");
        }
    }


    /**
     * Applies flags to the arguments and assign the default value
     * to result path if such argument is missing
     *
     * @param args
     * @return
     */
    private static void argSort(String[] args) {
//        String blockDictPath = null;
//        String genomePath = null;
//        String resultDir = null;
//        String extraParam = null;
//        String[] output = new String[4];
        try {
            for (int i = 0; i < args.length; i += 2) {
                String arg = args[i + 1];
                switch (args[i]) {
                    case SNPDICT:
                        snpDictPath = arg;
                        break;
                    case BLOCK:
                        blocksPath = arg;
                        break;
                    case RANGE:
                        SNPsearcher.range = Integer.parseInt(arg);
                        break;
                    case DELTA:
                        deltaFreq = Float.parseFloat(arg);
                        break;
                    default:
                        throw new Exception();
                }
            }
        } catch (Exception exp) {
            System.out.println("Invalid argument flag");
        }
//        // Default path for result file
//        if (resultDir == null) {
//            resultDir = ".";
//        }
    }
}

