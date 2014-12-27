import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import com.google.common.collect.*;

/**
* Static class containing tools to cyclize a linear Catalyst. 
*/
public final class MonteCarlo
{
    /** Desired bond length for forming bond */
    public static final double BOND_LENGTH = 1.35;

    /** Number of iterations; we always use this many. */
    public static final int ITERATIONS = 200;

    /** Initial temperature.  Final temp is zero. */
    public static final double KT = 2.0;

    /** Do not instantiate. */
    private MonteCarlo(){throw new IllegalArgumentException("Do not create instance of Monte Carlo!");}

    /** A molecule cyclizer.  Given a molecule and a list of rotatable bonds, the
    * method will make mutations to cyclize the endpoints to a normal bond length.
    * The return will be a new molecule, but with the same indices.
    * @param m the molecule to run this algorithm on
    * @param rotatableBonds the torsions that we can mutate
    * @param leftIndex the left side of the forming bond
    * @param rightIndex
    * @return a cyclized molecule; the new bond is not formed!
    */
    public static Molecule cyclize(Molecule m, List<IndexTorsion> rotatableBonds, int leftIndex, int rightIndex)
    {  
        double temperature = KT;

        System.out.println("Beginning Monte Carlo cyclization on " + m.name + ":\n:");

        for ( int i = 0 ; i < ITERATIONS ; i++ )
        {
            // this energy function includes the opls energy of the
            // linear fragment only, plus a triangular potential on the
            // terminal atoms.
            double oldEnergy = m.getOPLSenergy() + 100 * (Vector3D.distance(m.getAtom(leftIndex).position, m.getAtom(rightIndex).position) - BOND_LENGTH);
            Molecule testMolecule = mutate(m, rotatableBonds, temperature);
            double newEnergy = testMolecule.getOPLSenergy() + 100 * (Vector3D.distance(testMolecule.getAtom(leftIndex).position, testMolecule.getAtom(rightIndex).position) - BOND_LENGTH);
            if ( decider(newEnergy-oldEnergy, temperature) )
            {
                m = testMolecule;
                System.out.println("(" + i + ")Old Energy: " + oldEnergy + "\nNew Energy: " + newEnergy + "\n");
                System.out.println("{Accepted!}\n");
            }
           
            temperature = temperature - KT/ITERATIONS;
            
            /* for testing
            if ( i%100 == 0 )
                {
                    MOL2InputFile mol = new MOL2InputFile(m);
                    mol.write("opt_"+i/100+".mol2");
                }*/
        }
        
        return m;
    }

    /**
     * Takes an angle and restricts it to the range [-180.0, 180.0] using the modulus.
     * @param d an angle
     * @return the restricted angle
     */
    public static double angleModulus(double d) {
        double m = d % 360.0d;
        if (m < -180.0) m = m + 360.0;
        if (m >  180.0) m = m - 360.0;
        return m;
    }

    /** Decides whether to mutate.  Given a deltaE difference in energy and 
    * given temperature, it will return a Boltzmann distribution of true and 
    * false for unfavorable changes, true for favorable changes.
    */
    private static boolean decider(double deltaE, double temperature)
    {
        if ( deltaE < 0 )
            return true;
        else
        {
            if (ThreadLocalRandom.current().nextDouble() < Math.exp(-deltaE/temperature))
                return true;
            else
                return false;
        }
    }

    /** Mutation engine.  For each IndexTorsion, makes a random mutation of up to twenty degrees initially, but this amount gets smaller as the temperature decreases. */
    private static Molecule mutate(Molecule m, List<IndexTorsion> rotatableBonds, double temperature)
    {
        for ( IndexTorsion i : rotatableBonds )
        {  
            m = m.setDihedral(i,angleModulus(i.getDihedralAngle(m) - temperature/KT*(10 - (20 * ThreadLocalRandom.current().nextDouble()))));
        }
        return m;
    }
}

