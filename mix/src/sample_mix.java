import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class sample_mix {

    // Tissue files path
    private static ArrayList<String> matrices = new ArrayList<>();

    // Set of CpG sites that are not missing and common for all given tissues
    private static HashSet<String> goodSites = new HashSet<>();

    // Result dictionary with key - CpG site name, value - average frequency for each tissue for given CpG site
    private static SortedMap<String, ArrayList<String>> result;

    public static void main(String[] args) {
        // TODO: to write arg checker, that receives n samples with weight for each one, w1 + ... + wn = 1
        result = new TreeMap<>();
        matrices.addAll(Arrays.asList(args));
        filterSites();
        float[] coefficients = {0.47f, 0.05f, 0.47f};
        toMix(coefficients);

    }

    /**
     * Creates empty dictionary with CpG names
     */
    private static void prepareResult() {
        filterSites();
        ArrayList<String> resultLine = new ArrayList<>();
        for (String path : matrices) {
            String[] split = path.split("/");
            resultLine.add(split[split.length - 1].split("\\.")[0]);
        }
        result.put("cg", resultLine);
        for (String name : goodSites) {
            ArrayList<String> emptyList = new ArrayList<>();
            emptyList.add("0");
            result.put(name, emptyList);
        }
    }

    private static void toMix(float[] coefficients) {
        String line;
        String[] splitLine;
        int index = 0;
        try {
            prepareResult();
//            for (String tissue : matrices.keySet()) {
            for (String path : matrices) {
                BufferedReader tissue = getBuffer(path);
                line = tissue.readLine();
                line = tissue.readLine();
                while (line != null){
                    splitLine = line.split(",");
                    if (goodSites.contains(splitLine[0]) && !splitLine[1].equals("NA")) {
                        float newVal = Float.parseFloat(splitLine[1]) * coefficients[index];
                        result.get(splitLine[0]).set(0, Float.toString(Float.parseFloat(result.get(splitLine[0]).get(0)) + newVal));
                    }
                    line = tissue.readLine();
                }
                index ++;
            }
            saveMatrix();
        } catch (java.io.IOException exp) {
            System.out.println("Wrong file (sample_mix.toMix()");
        }
    }

    /**
     * Keeps only CpG sites that are common for all given tissues
     */
    private static void filterSites() {
        try {
            ArrayList<HashSet<String>> orgSets = new ArrayList<>();
            for (String path : matrices) {
                HashSet<String> preFilteredSet = new HashSet<>();
                BufferedReader matrix = getBuffer(path);
                String line = matrix.readLine();
                line = matrix.readLine();
                String[] splitLine;
                String name;
                while (line != null) {
                    splitLine = line.split(",");
                    name = splitLine[0];
                    preFilteredSet.add(name);
                    line = matrix.readLine();
                }
                orgSets.add(preFilteredSet);
            }
            if (orgSets.size() > 1) {
                HashSet<String> lastSet = orgSets.get(0);
                HashSet<String> big;
                HashSet<String> small;
                for (int i = 1; i < orgSets.size(); i++) {
                    if (orgSets.get(i).size() > lastSet.size()) {
                        big = orgSets.get(i);
                        small = lastSet;
                    } else {
                        big = lastSet;
                        small = orgSets.get(i);
                    }
                    big.retainAll(small);
                    lastSet = big;
                }
                orgSets.clear();
                orgSets.add(lastSet);
            }
            goodSites = orgSets.get(0);

        } catch (java.io.IOException exp) {
            System.out.println("Wrong file (sample_mix.filterSites())");
        }
    }


    /**
     * Creates new buffer reader for given file
     * @param path path to the file
     * @return new buffer reader for given file
     */
    private static BufferedReader getBuffer(String path) {
        try {
            InputStream fileStream = new FileInputStream(path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader bdecoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
            return new BufferedReader(bdecoder);
        } catch (java.io.IOException exp) {
            System.out.println("Wrong file (sample_mix.getBuffer())");
        }
        return null;
    }

    /**
     * Saves result matrix into the .csv file in project directory
     */
    private static void saveMatrix() {
        try {
//            Path resultPath = Paths.get(System.getProperty("user.dir"), "/result_file.csv");
            Path resultPath = Paths.get("/cs/cbio/ilia/data/project_data/tissue_samples/mixtures/skin_sperm_urine/Skin_47_Sperm_5_Urine_47.csv");
            PrintWriter out = new PrintWriter(resultPath.toString());
            for (String name : result.keySet()) {
                ArrayList <String> line = result.get(name);
                out.print(name + ",");
                for (int i = 0; i < line.size() - 1; i++) {
                    out.print(line.get(i) + ",");
                }
                out.print(line.get(line.size() - 1) + "\n");
            }
            out.close();
        } catch (java.io.IOException exp) {
            System.out.println("Matrix saving problem");
        }
    }
}
