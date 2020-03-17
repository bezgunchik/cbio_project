import java.util.Comparator;

public class GenomePlaceComparator implements Comparator<String[]> {


    @Override
    public int compare(String[] t1, String[] t2) {
        // Currently doesn't support chromosomes X and Y
        if (t1[0].equals(t2[0])) {
            return Integer.compare(Integer.parseInt(t1[1]), Integer.parseInt(t2[1]));
        } else {
            return Integer.compare(Integer.parseInt(t1[0].substring(3)), Integer.parseInt(t2[0].substring(3)));
        }
    }
}
