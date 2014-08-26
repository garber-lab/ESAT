package umms.core.coordinatesystem;

public class PermutationNotFoundException extends RuntimeException {
    public PermutationNotFoundException(int n)
    {
       super("Unable to find a permutation match after " + n + " tries.");
    }
}
