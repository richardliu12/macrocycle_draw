import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
* Handles the analysis of a .mae output.  This output file should contain
* multiple structures from a minimization of macrocycle conformations.
* The class will read in these geometries and energies, then produce
* .gjf files of the conformations.  Also, a list of the structures by
* energy and by inter(thio)urea distances will be produced.
*/
public abstract class Analysis
{
    public static void main(String[] args)
    {
        String filename = args[0];
        // read the .mae from the appropriate folder
        MAEFile m = new MAEFile(Settings.WORKING_DIRECTORY + "output/" + filename + "/" + filename + "-csearch-min.mae");

        // deserialize catalyst
        Catalyst c = null;
        try
        {
            FileInputStream fileIn = new FileInputStream(Settings.WORKING_DIRECTORY + 
                "output/" + filename + "/catalyst.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            c = (Catalyst) in.readObject();
            in.close();
            fileIn.close();
        }
        catch(IOException i)
        {
            i.printStackTrace();
        }
        catch(ClassNotFoundException e)
        {
            System.out.println("Unable to find catalyst.ser!");
        }

        int counter = 0;
        for ( List<Vector3D> geom : m.geometries )
            {
                // map the new conformation to the catalyst and print out gjfs
                Map<Atom,Atom> atomMap = new HashMap<>();
                if ( c.contents.size() != geom.size())
                    throw new IllegalArgumentException("Sizes of molecules do not match!");
                for ( int i = 0; i < c.contents.size(); i++ )
                    atomMap.put(c.getAtom(i+1), c.getAtom(i+1).moveAtom(geom.get(i)));
                c = c.moveAtoms(atomMap);

                String moleculeName = c.name + String.format("%05d", ++counter);
                GaussianInputFile file = new GaussianInputFile(c, c.name, 
                    "#p opt b3lyp m06-2x geom=connect scrf=(solvent=benzene) freq");
                file.write(Settings.WORKING_DIRECTORY + "output/" + filename + "/" + moleculeName + ".gjf");
            }
    }
}
