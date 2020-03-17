import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;

//zcat Skin.csv.gz | grep -E "cg[0-9]+(,NA){2,}"

public class matrix_average {

    // Tissue files path
    private ArrayList<String> matrices;

    // Set of CpG sites that are totally missing in at least one of the tissues
    private HashSet<String> missingSites;

    // Set of CpG sites that are not missing and common for all given tissues
    private HashSet<String> goodSites;

    // Result dictionary with key - CpG site name, value - average frequency for each tissue for given CpG site
    private SortedMap<String, ArrayList<String>> result;


    /**
     * Constructor of matrix_average class.
     * Saves tissue files paths and creates a file with CpG that are totally missing in at least one of the tissues
     * @param filePaths tissue files paths
     */
    matrix_average(String[] filePaths) {
        matrices = new ArrayList<>();
        missingSites = new HashSet<>();
        goodSites = new HashSet<>();
        result = new TreeMap<>();
        try {
            Path tempFilePath = Paths.get(System.getProperty("user.dir"), "/temp_file" + ".csv");
            Runtime.getRuntime().exec("touch " + tempFilePath.toString());
            for (String path : filePaths) {
                matrices.add(path);
                String[] numCmd = {"sh",
                        "-c", "zcat " + path + " | awk '{print split($0,a,\",\"); exit}'"};
                int colNum = Integer.parseInt(new BufferedReader(new InputStreamReader(Runtime.getRuntime().
                        exec(numCmd).getInputStream())).readLine()) - 1;
                String[] command = {"sh", "-c", "zcat " + path + " | grep -E \"cg[0-9]+(,NA){" +
                        colNum + "}\" | awk '{split($0,a,\",\"); print a[1]}' >> " + tempFilePath.toString()};
                Runtime.getRuntime().exec(command);
            }
            BufferedReader fileMissing = new BufferedReader(new FileReader(tempFilePath.toString()));
            String line = fileMissing.readLine();
            while (line != null) {
                missingSites.add(line.split(",")[0]);
                line = fileMissing.readLine();
            }
//            Runtime.getRuntime().exec("rm -f " + tempFilePath.toString());
            fileMissing.close();
        } catch (java.io.IOException exp) {
//            System.out.println("Wrong file (matrix_average constructor)");
            System.out.println(exp.getMessage());
        }
    }


    /**
     * Keeps only CpG sites that are common for all given tissues
     */
    private void filterSites() {
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
                    if (!missingSites.contains(name)) {
                        preFilteredSet.add(name);
                    }
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
            System.out.println("Wrong file (matrix_average.getAverage())");
        }
    }


    /**
     * Creates empty dictionary with CpG names
     */
    private void prepareResult() {
        filterSites();
        ArrayList<String> resultLine = new ArrayList<>();
        for (String path : matrices) {
            String[] split = path.split("/");
            resultLine.add(split[split.length - 1].split("\\.")[0]);
        }
        result.put("cg", resultLine);
        for (String name : goodSites) {
            result.put(name, new ArrayList<>());
        }
    }


    /**
     * Calculates average frequency for each tissue for each CpG site
     */
    void getAverage() {
        String line;
        String[] splitLine;
        try {
            prepareResult();
//            for (String tissue : matrices.keySet()) {
            for (String path : matrices) {
                BufferedReader tissue = getBuffer(path);
                line = tissue.readLine();
                line = tissue.readLine();
                while (line != null){
                    splitLine = line.split(",");
                    float counter = 0;
                    float average = 0;
                    if (goodSites.contains(splitLine[0])) {
                        for (String value : Arrays.copyOfRange(splitLine, 1, splitLine.length)) {
                            if (!value.equals("NA")) {
                                counter ++;
                                average += Float.parseFloat(value);
                            }
                        }
                        result.get(splitLine[0]).add(Float.toString(average/counter));
                    }
                    line = tissue.readLine();
                }
            }
        saveMatrix();
        } catch (java.io.IOException exp) {
            System.out.println("Wrong file (matrix_average.getAverage()");
        }
    }


    /**
     * Creates new buffer reader for given file
     * @param path path to the file
     * @return new buffer reader for given file
     */
    private BufferedReader getBuffer(String path) {
        try {
            InputStream fileStream = new FileInputStream(path);
            InputStream gzipStream = new GZIPInputStream(fileStream);
            Reader bdecoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
            return new BufferedReader(bdecoder);
        } catch (java.io.IOException exp) {
            System.out.println("Wrong file (matrix_average.getAverage()");
        }
        return null;
    }


    /**
     * Saves result matrix into the .csv file in project directory
     */
    private void saveMatrix() {
        try {
//            Path resultPath = Paths.get(System.getProperty("user.dir"), "/result_file.csv");
            Path resultPath = Paths.get("/cs/cbio/ilia/data/project_data/tissue_samples/matrices/avrg_matrix_minus_last_col2.csv");
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

///cs/cbio/ilia/data/project_data/tissue_samples/matrices/Urine_minus_last_col.csv.gz /cs/cbio/ilia/data/project_data/tissue_samples/matrices/Sperm_minus_last_col.csv.gz /cs/cbio/ilia/data/project_data/tissue_samples/matrices/Skin_minus_last_col.csv.gz /cs/cbio/ilia/data/project_data/tissue_samples/matrices/Saliva_minus_last_col.csv.gz /cs/cbio/ilia/data/project_data/tissue_samples/matrices/Whole_blood_minus_last_col.csv.gz
// TODO: Сделать несколько атласов, в которых отсутствует по одному экземпляру и сохранить этот экземпляр отдельно
// TODO: (в парах атлас-экземпляр)



//TODO: понять, где добавляются лишние CpG