import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * A class that is used to parse the output from a Tinker call to analyze.
 * It creates an energy breakdown of the total energy of a molecule into residue
 * components by treating the energy of a component as a self residue and an
 * interaction term.
 **/
public class GaussianMoleculeOutputFile extends OutputFileFormat
{

    /** The total energy of the molecule **/
    public final double totalEnergy;
    /** A unique identifying number **/
    public final int jobID;
    
    public String toString()
    {
        return "Total Energy:  " + totalEnergy;
    }

    /** fileName - the name of the file containing output from a Tinker call to analyze 
	molecule = the molecule of the corresponding analyze file whose energy will be divided up by residue
    **/
    public GaussianMoleculeOutputFile(String fileName, int id)
    {
        super(fileName);
        jobID = id;
        // sum of all individual energies used to check if method works
        Double totalEnergy = 0.0;


        for (List<String> line : fileContents) 
        {
            try
            {
                if (line.get(0).equals("SCF") && line.get(1).equals("Done:"))
                {
                    totalEnergy = Double.parseDouble(line.get(4));
                    break;
                }
            }
            catch (Exception e)
            {
                //System.out.println("filename: " + fileName);
                //System.out.println("current line:");
                //System.out.print(line);
                throw e;
            }
        }

        // convert to kCal/mol
        this.totalEnergy = totalEnergy*627.509469;
    } // end of constructor

    private int getInt(String s)
    {
        String newString = s.replaceAll("[^\\d]","");
        return Integer.parseInt(newString);
    }

    /** for testing */
    public static void main(String[] args) 
    {
        GaussianMoleculeOutputFile test = new GaussianMoleculeOutputFile("test.out", 1);
        System.out.println(test.toString());
    }
} // end of class GaussianMoleculeOutputFile
