import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
* Static class containing tools to cyclize a linear Catalyst. 
*/
public final class MonteCarlo
{
    /** Desired bond length */
    public static final double BOND_LENGTH = 1.35;

    /** Number of iterations */
    public static final int ITERATIONS = 10000;

    /** Initial temperature */
    public static final double KT = 1;

    /** Do not instantiate. */
    private MonteCarlo(){}

    /** A molecule cyclizer.  Given a molecule and a list of rotatable bonds, the
    * method will make mutations to cyclize the endpoints to a normal bond length.
    * The return will be a new molecule, but with the same indices.
    */
    public static Molecule cyclize(Molecule m, List<IndexTorsion> rotatableBonds, int leftIndex, int rightIndex)
    {  
        double temperature = KT;

        for ( int i = 0 ; i < ITERATIONS ; i++ )
        {
            double oldEnergy = m.getOPLSenergy() + 1000 * (Vector3D.distance(m.getAtom(leftIndex).position, m.getAtom(rightIndex).position) - BOND_LENGTH);
            Molecule testMolecule = mutate(m, rotatableBonds);
            double newEnergy = testMolecule.getOPLSenergy() + 1000 * (Vector3D.distance(testMolecule.getAtom(leftIndex).position, testMolecule.getAtom(rightIndex).position) - BOND_LENGTH);
            if ( decider(newEnergy-oldEnergy, temperature) )
            {
                m = testMolecule;
                System.out.println("{Accepted!}");
            }

            System.out.println("Old Energy: " + oldEnergy + "\n\nNew Energy: " + newEnergy);
            
            temperature = temperature*(1 - 1/ITERATIONS);
        }
        
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
            if (Math.random() < Math.exp(-deltaE/temperature))
                return true;
            else
                return false;
        }
    }

    /** Mutation engine.  For each IndexTorsion, makes a random mutation of up to ten degrees. */
    private static Molecule mutate(Molecule m, List<IndexTorsion> rotatableBonds)
    {
        for ( IndexTorsion i : rotatableBonds )
        {  
            m = m.setDihedral(i, i.getDihedralAngle(m) - 5 + (10 * Math.random()));
        }
        return m;
    }
}

