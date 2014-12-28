import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * The main program.  Draws macrocycles from input files and given template.
 * These files should be in ./input.  These cyclized catalysts will be
 * minimized and a .mae file will be put in ./mae, as well as a .mol2 file.
 * Next, a conformational search file is created for each catalyst, as well
 * as a serialized Catalyst file.  Once executed, the results can be found
 * ./output.  This output should be analyzed using a separate script.
 */
public class Macrocycle_Draw
{
    public static void main(String[] args)
    {
        // read from file "template" in the input directory to get
        // the catalyst template that we want to use
        File file = new File(Settings.INPUT_DIRECTORY + "template");
        String firstLine = null; // we will only read the first line
        boolean C2 = false; // are we in C2 mode?
        List<FragmentType> template = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(Settings.INPUT_DIRECTORY + "template")))
            {
                String line = reader.readLine();
                String[] fields = line.split("@");
                for ( String s : fields )
                    {
                        switch(s.toLowerCase())
                            {
                                case "c2": C2 = true;
                                    break;
                                case "urea": template.add(FragmentType.UREA);
                                    break;
                                case "thiourea": template.add(FragmentType.THIOUREA);
                                    break;
                                case "linker_1": template.add(FragmentType.LINKER_1);
                                    break;
                                case "linker_2": template.add(FragmentType.LINKER_2);
                                    break;
                                case "linker_3": template.add(FragmentType.LINKER_3);
                                    break;
                                case "linker_4": template.add(FragmentType.LINKER_4);
                                    break;
                                case "": break;
                                default: throw new IllegalArgumentException("Error in template file, cannot read: " + s);
                           }
                    }
            }
        catch (Exception e)
            {
               throw new IllegalArgumentException(e.getMessage());
            }

        System.out.println("Reading template... \n" + template);

        // Create the catalyst list.  These have been cyclized, but 
        // not yet minimized.
        List<Catalyst> catalysts = new ArrayList<>();

        if (C2) 
            catalysts = FragmentLibrary.createC2Catalysts(template);
        else
            catalysts = FragmentLibrary.createCatalysts(template);

        // Creates mol2 files for each of the catalysts
        for ( Catalyst c : catalysts )
            {
                MOL2InputFile m = new MOL2InputFile(c);
                m.write(c.name + ".mol2");
            }

        // Minimization
        for ( Catalyst c : catalysts )
            {
                System.out.println("Minimizing " + c.name + "...");
                try
                    {
                        Runtime runtime = Runtime.getRuntime();
                        String runString = Settings.WORKING_DIRECTORY + "minimization.sh " + c.name;
                        Process process = runtime.exec(runString);
                        int exitValue = process.waitFor();
                    }
                catch (Exception e)
                    {
                        System.out.println("Error during minimization job:");
                        e.printStackTrace();
                    }
                System.out.println("DONE!");
            }
    
        // Writes a .com file for a catalyst
        for ( Catalyst c : catalysts )
            {   
                System.out.println("Writing conformational search file for " +
                    c.name + "...");
                COMInputFile cfile = new COMInputFile(c);
                cfile.write(Settings.WORKING_DIRECTORY + c.name + ".com");
                System.out.println("DONE!\n");
            }
    }
}
