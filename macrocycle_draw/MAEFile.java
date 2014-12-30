import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * Parses Macromodel geometry files.  It's expected a full block and partial blocks will be present.
 * I'm not bothering to read the atoms or any other data.
 */
public class MAEFile extends OutputFileFormat
{
    /** The energies of the geometries, converted from kJ to kcal. */
    public final List<Double> energies;

    /** The RMS gradients of the minimized structures. */
    public final List<Double> gradients;

    /** The geometries in the file. */
    public final List<List<Vector3D>> geometries;

    /**
     * Reads Macromodel energies, gradients, and geometries from the specified filename.
     */
    public MAEFile(String filename)
    {
        super(filename);
        
        // make some temporary fields
        List<Double> tempEnergies = new ArrayList<>();
        List<Double> tempGradients = new ArrayList<>();
        List<List<Vector3D>> tempGeometries = new ArrayList<>();

        // get blocks
        StringBuilder fullBlock = new StringBuilder();
        StringBuilder partialBlock = new StringBuilder();

        // keeps track of which section of the file we are in
        boolean newBlock = false;
        boolean pastHeader = false;
        boolean readPartialBlocks = false;
        List<List<String>> currentBlock = new ArrayList<>();

        for (int i=0; i < fileContents.size(); i++)
            {
                // get current line
                List<String> fields = fileContents.get(i);
                //System.out.println(fields);

                // figure out what part of the file we're in
                if ( !pastHeader )
                    {
                        if ( fields.get(0).indexOf("f_m_ct") > -1 )
                            pastHeader = true;
                        else
                            continue;
                    }

                if ( fields.get(0).indexOf("p_m_ct") > -1 || i == fileContents.size() - 1 )
                    {
                        if ( !readPartialBlocks )
                            {
                                parseFullBlock(currentBlock, tempEnergies, tempGradients, tempGeometries);
                                readPartialBlocks = true;
                            }
                        else
                            parsePartialBlock(currentBlock, tempEnergies, tempGradients, tempGeometries);
                        currentBlock = new ArrayList<>();
                    }
                else
                    currentBlock.add(fields);
            }

        // set fields
        this.energies = ImmutableList.copyOf(tempEnergies);
        this.gradients = ImmutableList.copyOf(tempGradients);
        this.geometries = ImmutableList.copyOf(tempGeometries);
    }

    /**
     * Debugging method that prints out a block of text.
     */
    public static void printBlock(List<List<String>> block)
    {
        System.out.println("================");
        for (List<String> line : block)
            {
                String newString = "";
                for (String s : line)
                    newString += s + " ";
                System.out.printf("%-50s\n", newString);
            }
        System.out.println("================");
    }

    /**
     * Parses a "f_m_ct" block and updates the lists that are passed in.
     * @param block the block of text to parse
     * @param energies where to update the energies to
     * @param gradient where to update the gradients to
     * @param geometries where to update the geometries to
     */
    public static void parseFullBlock(List<List<String>> block, List<Double> energies, List<Double> gradients, List<List<Vector3D>> geometries)
    {
        Double energy = Double.valueOf(block.get(14).get(0)) / 4.184;
        Double gradient = Double.valueOf(block.get(23).get(0));
        List<Vector3D> geometry = new ArrayList<>();
        int section = 0;
        for ( int i=0; i < block.size(); i++ )
            {
                List<String> fields = block.get(i);
                if ( fields.get(0).indexOf(":::") > -1 )
                    {
                        section++;
                        continue;
                    }

                if ( section == 4 )
                    {
                        double x = Double.parseDouble(fields.get(2));
                        double y = Double.parseDouble(fields.get(3));
                        double z = Double.parseDouble(fields.get(4));
                        geometry.add(new Vector3D(x,y,z));
                    }
            }
        energies.add(energy);
        gradients.add(gradient);
        geometries.add(ImmutableList.copyOf(geometry));
   }
    
    /**
     * Parses a "p_m_ct" block and updates the lists that are passed in.
     * @param block the block of text to parse
     * @param energies where to update the energies to
     * @param gradient where to update the gradients to
     * @param geometries where to update the geometries to
     */
    public static void parsePartialBlock(List<List<String>> block, List<Double> energies, List<Double> gradients, List<List<Vector3D>> geometries)
    {
        Double energy = Double.valueOf(block.get(13).get(0)) / 4.184;
        Double gradient = Double.valueOf(block.get(22).get(0));
        List<Vector3D> geometry = new ArrayList<>();
        int section = 0;
        for ( int i=0; i < block.size(); i++ )
            {
                List<String> fields = block.get(i);
                if ( fields.get(0).indexOf(":::") > -1 )
                    {
                        section++;
                        continue;
                    }

                if ( section == 4 )
                    {
                        double x = Double.parseDouble(fields.get(1));
                        double y = Double.parseDouble(fields.get(2));
                        double z = Double.parseDouble(fields.get(3));
                        geometry.add(new Vector3D(x,y,z));
                    }
            }
        energies.add(energy);
        gradients.add(gradient);
        geometries.add(ImmutableList.copyOf(geometry));
    }

    @Override
    public String toString()
    {
        return String.format("MAEFile (%d structures)", geometries.size());
    }

    /** for testing */
    public static void main(String[] args) 
    {
        MAEFile file = new MAEFile("test.mae");
        List<List<Vector3D>> geometries = file.geometries;
        for (Double energy : file.energies)
            System.out.println(energy);
    }
}
