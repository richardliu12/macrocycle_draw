import java.io.*;
import java.util.*;
import org.jgrapht.alg.*;
import org.jgrapht.graph.*;

/**
 * This class manages output to a .com file for MacroModel conformational search.
 * The class itself only stores the associated catalyst, but
 * the write method will write the parameters necessary to 
 * run a conformational search.  See the main program class
 * for organization of input and output files.
 */
public class COMInputFile // does not extend InputFileFormat
{
    /** The associated Catalyst.  Must have already been cyclized
     * with .mae file in the ./mae folder.
     */
    public final Catalyst catalyst;

    /**
     * Constructs the .com file from a Catalyst.  Actually, this
     * does not do anything until you write the file to disk.
     */
    public COMInputFile(Catalyst catalyst)
    {
        this.catalyst = catalyst;
    }

    public void write(String filename)
    {
        // First copy the header from the template file
        try
        {
            Runtime runtime = Runtime.getRuntime();
            String runString = Settings.WORKING_DIRECTORY + "make_csearch.sh "
                + catalyst.name;
            Process process = runtime.exec(runString);
            int exitValue = process.waitFor();
        }
        catch (Exception e)
        {
            System.out.println("Error during creation of conf. search: ");
            e.printStackTrace();
        }

        //
        // Find comparison atoms, and add to file
        //
        List<Integer> heavyAtomIndices = new ArrayList<>();
        for ( Fragment f : catalyst.fragmentList )
            for ( Atom a : f.contents )
                if ( a.element != Element.HYDROGEN && a.element != Element.DUMMY )
                    heavyAtomIndices.add(catalyst.getAtomNumber(a));

            // Make the heavyAtom list into a multiple of 4 for neatness
        while ( heavyAtomIndices.size()%4 != 0 )
            heavyAtomIndices.add(0);
        
        for ( int i = 0; i < heavyAtomIndices.size(); i+=4 )
            {
                String writeString = toCOMLine("COMP", heavyAtomIndices.get(i),
                    heavyAtomIndices.get(i+1), heavyAtomIndices.get(i+2),
                    heavyAtomIndices.get(i+3), 0, 0, 0, 0);
                InputFileFormat.appendStringToDisk("\n" + writeString, 
                    Settings.WORKING_DIRECTORY + "/mae/"
                    + catalyst.name + ".com");
            }

        //
        // Find chiral atoms, and add to file
        //
        List<Integer> chiralAtomIndices = new ArrayList<>();
        for ( Fragment f : catalyst.fragmentList )
            for ( Atom a : f.chiralAtoms )
                chiralAtomIndices.add(catalyst.getAtomNumber(a));

        for ( int i = 0; i < chiralAtomIndices.size(); i++)
            InputFileFormat.appendStringToDisk("\n" + toCOMLine("CHIG",
                chiralAtomIndices.get(i), 0, 0, 0, 0, 0, 0, 0), Settings.WORKING_DIRECTORY
                + "/mae/" + catalyst.name + ".com");

        //
        // Find torsions, and add to file
        //
        List<IndexTorsion> torsions = catalyst.getLinearTorsions();
        for ( IndexTorsion i : torsions )
            InputFileFormat.appendStringToDisk("\n" + toCOMLine("TORS",
                i.index2, i.index3, 0, 0, 0, 180, 0, 0), Settings.WORKING_DIRECTORY
                + "/mae/" + catalyst.name + ".com");

        //
        // Add the ring closure
        //
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> tempConnectivity = catalyst.connectivity;
        Atom leftConnect = catalyst.fragmentList.get(0).leftConnect;
        Atom rightConnect = catalyst.fragmentList.get(catalyst.fragmentList.size()-1).rightConnect;
        tempConnectivity.removeEdge(leftConnect,rightConnect);
        List<DefaultWeightedEdge> shortestPath = DijkstraShortestPath.findPathBetween(tempConnectivity, leftConnect, rightConnect);

        Atom leftConnectNext = null, rightConnectNext = null;

            // search for the ring atoms that are connected to leftConnect and rightConnect
        for ( DefaultWeightedEdge e : shortestPath )
            {
                if ( tempConnectivity.getEdgeSource(e) == leftConnect )
                    leftConnectNext = tempConnectivity.getEdgeTarget(e);
                if ( tempConnectivity.getEdgeTarget(e) == leftConnect )
                    leftConnectNext = tempConnectivity.getEdgeSource(e);
                if ( tempConnectivity.getEdgeSource(e) == rightConnect )
                    rightConnectNext = tempConnectivity.getEdgeTarget(e);
                if ( tempConnectivity.getEdgeTarget(e) == rightConnect )
                    rightConnectNext = tempConnectivity.getEdgeSource(e);
            }
        InputFileFormat.appendStringToDisk("\n" + toCOMLine("RCA4", catalyst.getAtomNumber(leftConnectNext),
            catalyst.getAtomNumber(leftConnect), catalyst.getAtomNumber(rightConnect),
            catalyst.getAtomNumber(rightConnectNext), 0.5, 2.5, 0, 0),
            Settings.WORKING_DIRECTORY + "/mae/" + catalyst.name + ".com");
                  
        //
        // Add the footer
        //
        InputFileFormat.appendStringToDisk("\n" + toCOMLine("CONV", 2, 0, 0, 0,
            0.5, 0, 0, 0), Settings.WORKING_DIRECTORY + "/mae/" + catalyst.name
            + ".com");
        InputFileFormat.appendStringToDisk("\n" + toCOMLine("MINI", 1, 0, 1000, 0,
            0, 0, 0, 0), Settings.WORKING_DIRECTORY + "/mae/" + catalyst.name + ".com");
    
        //
        // Serialize the Catalyst file
        //
        try
        {
            new File(Settings.WORKING_DIRECTORY + "/output/" + catalyst.name).mkdir();
            FileOutputStream fileOut = new FileOutputStream(Settings.WORKING_DIRECTORY + 
                "/output/" + catalyst.name + "/catalyst.ser");
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(this.catalyst);
            out.close();
            fileOut.close();
        }
        catch(IOException i)
        {
            i.printStackTrace();
        }
    }

    /** Formats a line for .com file.
     * @param command
     * @param i1
     * @param i2
     * @param i3
     * @param i4
     * @param d1
     * @param d2
     * @param d3
     * @param d4
     * @return the line
     */
    public static String toCOMLine(String command, int i1, int i2, int i3, int i4, double d1, double d2, double d3, double d4)
    {
        String returnString = String.format("%4s%8d%7d%7d%7d%11.4f%11.4f%11.4f%11.4f", 
                    command.toUpperCase(), i1, i2, i3, i4, d1, d2, d3, d4);
        return (" " + returnString); // initial space is crucial
    }
}
