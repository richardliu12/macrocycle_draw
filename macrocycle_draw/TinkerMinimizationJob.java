import java.util.*;
import java.io.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;


/**
 * This class represents a TINKER minimization job.
 * It is assumed that TINKER is in the path.  On Cygwin, check that the
 * executable files end in .exe.
 */
public class TinkerMinimizationJob implements WorkUnit, Serializable, Immutable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** help us choose a unique filename */
    public static final AtomicInteger index = new AtomicInteger();

    /** the xyz file that tinker will read from */
    public final TinkerXYZInputFile tinkerXYZInputFile;

    /** the key file that tinker will read from */
    public final TinkerKeyFile tinkerKeyFile;

    /**
     * Creates a job for minimizing some molecule.
     * @param molecule the molecule whose geometry is to be minimized
     * @param extraKeywords any extra keywords that are desired (newlines required)
     */
    public TinkerMinimizationJob(Molecule molecule, String extraKeywords)
    {
        if ( molecule == null )
            throw new NullPointerException("molecule cannot be null when constructing a tinker minimization job");
       	tinkerXYZInputFile = new TinkerXYZInputFile(molecule);
        String keywords = Settings.TINKER_MINIMIZATION_STANDARD_KEYWORDS + extraKeywords;
	    tinkerKeyFile = new TinkerKeyFile(keywords);
    }
    
    /*
    /**
     * Creates a job for minimizing a peptide pose based on a template peptide
     * @param templatePeptide the template peptide in the non-close contact form
     * @param pose the positions of the atoms in a specific pose (either ground or transition state)
     * @param extraKeywords any extra keywords that are desired (newlines required)
     */
     /*
     public TinkerMinimizationJob(Peptide templatePeptide, List<Vector3D> pose, String extraKeywords)
     {
         if (templatePeptide.contents.size() != pose.size())
             throw new IllegalArgumentException("Size of template peptide's contents must match the pose size");
        
        //creates String to write
        String outputString = "";
        //write number of atoms and molecule name
        outputString = templatePeptide.contents.size() + " " + templatePeptide.name + "\n";

        //write atom list and connections
        for (int currentAtomNumber = 1; currentAtomNumber <= templatePeptide.contents.size(); currentAtomNumber++)
            {
                Atom currentAtom = templatePeptide.contents.get(currentAtomNumber);
                Set<DefaultWeightedEdge> bonds = templatePeptide.connectivity.edgesOf(currentAtom);
                outputString = outputString + String.format("%3d %2s %12.8f %12.8f %12.8f %6d", currentAtomNumber,
                                                            currentAtom.element.symbol,  pose.get(currentAtomNumber).getX(),
                                                            pose.get(currentAtomNumber).position.getY(), pose.get(currentAtomNumber).getZ(),
                                                            currentAtom.tinkerAtomType);

                for (DefaultWeightedEdge e : bonds)
                    {
                        int edgeTarget = contents.indexOf(templatePeptide.connectivity.getEdgeTarget(e)) + 1;
                        int edgeSource = contents.indexOf(templatePeptide.connectivity.getEdgeSource(e)) + 1;
                        int edgeTargetToWrite = edgeTarget;
                        //change to edgeSource if edgeTarget is simply the current atom (edges have no particular ordering)
                        if (edgeTarget == currentAtomNumber)
                            edgeTargetToWrite = edgeSource;
                        outputString = outputString+String.format("%6d",edgeTargetToWrite);
                    }
                outputString = outputString + "\n";
            }
       
        
        tinkerXYZInputFile = new TinkerXYZInputFile(outputString);
        String keywords = Settings.TINKER_MINIMIZATION_STANDARD_KEYWORDS + extraKeywords;
	    tinkerKeyFile = new TinkerKeyFile(keywords);
    }
    */

    /** will run with standard keywords only */
    public TinkerMinimizationJob(Molecule molecule)
    {
        this(molecule,"");
    }

    /**
     * Auto-selects a filename and runs the minimization calculation.
     * It is assumed that the TINKER binaries are in the path.
     * @return the result of the calculation
     */
    public TinkerMinimizationResult call()
    {
        // check whether the shared memory folder is there
        // guard against simultaneous copying
        synchronized (Settings.SHM_LOCK)
            {
                File testFile = new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY);
                if ( ! testFile.exists() )
                    {
                        // if the shared memory directory isn't there, copy it over
                        System.out.print("Shared memory folder deleted, so copying it again...");
                        String runString = Settings.WORKING_DIRECTORY + "copyshm.sh";
                        try
                            {
                                Process process = Runtime.getRuntime().exec(runString);
                                process.waitFor();
                            }
                        catch (Exception e)
                            {
                                e.printStackTrace();
                            }
                        System.out.println("done.");
                    }
            }

        // choose an appropriate base filename
        // try up to TINKER_MINIMIZATION_MAX_FILENAMES times to make a unique set of filenames
        String baseFilename = "";
        counting:
        for (int i=0; i < Settings.TINKER_MINIMIZATION_MAX_FILENAMES; i++)
            {
                // get a new ID number for this job
                int currentIndex = index.getAndIncrement();
                baseFilename = String.format("%s_tinker_minimization_job_%010d", Settings.HOSTNAME, currentIndex);

                // don't allow this choice of filenames if any files with this prefix already exist
                for ( File f : new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles() )
                    {
                        if ( f.getName().startsWith(baseFilename) )
                            {
                                baseFilename = "";
                                continue counting;
                            }
                    }
                break;
            }
        if ( baseFilename.length() == 0 )
            throw new IllegalArgumentException("Unable to set filename!");

        // reset counter if necessary
        if ( index.get() > Settings.TINKER_MINIMIZATION_MAX_FILENAMES )
            index.getAndSet(0);

	    // write input files to disk
        tinkerXYZInputFile.write( Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".xyz" );
	    tinkerKeyFile.write     ( Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".key" );

        // call minimize
        double elapsedTime = 0.0;
        int exitValue = -1;
        //boolean badGeometry = false;
        try
            {
                if ( Settings.PLATFORM == Settings.Platform.DOS )
                    {
                        ProcessBuilder builder = new ProcessBuilder("cmd", "/c", "run_tinker.bat", baseFilename);
                        builder.directory(new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY));
                        long startTime = System.currentTimeMillis();
                        Process process = builder.start();
                        exitValue = process.waitFor();
                        //System.out.println(exitValue);
                        long endTime = System.currentTimeMillis();
                        elapsedTime = (endTime - startTime) / 1000.0;
		            }
                else if ( Settings.PLATFORM == Settings.Platform.LINUX )
                    {
                        long startTime = System.currentTimeMillis();
                        String runString = Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + "run_tinker.sh " +
                                           Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + " " + baseFilename;
                        //System.out.println(runString);
                        Process process = Runtime.getRuntime().exec(runString);
                        
                        while ( true )
                            {
                                // check if the process has terminated
                                GeneralThreadService.wait(500);
                               
                                /*try
                                    {
                                        TinkerMinimizationLogFile tempOutput = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
                                        if ( tempOutput.stringRepresentation.indexOf("Incomplete Convergence due to BadIntpln") > -1 )
                                            badGeometry = true;
                                    }
                                catch (Exception e) {}*/

                                boolean done = true;
                                try { exitValue = process.exitValue(); }
                                catch (Exception e) { done = false; }
                                //if ( done || badGeometry )
                                if ( done )
                                    break;

                                long now = System.currentTimeMillis();
                                elapsedTime = (now - startTime)/1000.0;
                                //if ( elapsedTime > 45.0 || badGeometry )
                                if ( elapsedTime > 120.0 )
                                    {
                                        process.destroy();
                                        break;
                                    }
                            }
                        long endTime = System.currentTimeMillis();
                        elapsedTime = (endTime - startTime) / 1000.0;
                    }
            }
        catch (Exception e)
            {
                System.out.println("Error while running Tinker job:");
                System.out.println(baseFilename);
                e.printStackTrace();
            }

        // remove files on abnormal termination
        if ( exitValue != 0 )
            {
                try
                    {
                        File[] files = new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles();
                        for ( File f : files )
                            {
                                String filename = f.getName();
                                if ( f.getName().startsWith(baseFilename) )
                                    {
                                        //f.delete();
                                        //System.out.println(f.getName() + " deleted.");
                                    }
                            }
                    }
                catch (Exception e)
                    {
                        System.out.println("Error while trying to delete files:");
                        e.printStackTrace();
                    }
            }

        // check if the job completed correctly
        if ( exitValue == -1 )
            throw new IllegalArgumentException(baseFilename + " exceeded the allotted time");
        //else if ( badGeometry )
        //    throw new IllegalArgumentException(baseFilename + ": interpolation error");
        else if ( exitValue != 0 )
            {
                String tail = "";
                try
                    {
                        TinkerMinimizationLogFile errorOutput = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
                        String[] lines = errorOutput.stringRepresentation.split("\n");
                        int length = lines.length;
                        for (int i=Math.max(length-10,0); i < length; i++)
                            tail += lines[i] + "\n";
                    }
                catch (Exception e)
                    {
                    }
                throw new IllegalArgumentException("error code " + exitValue + " while runing tinker minimization job: " + baseFilename + "!\n" + tail);
            }

        // retrieve output XYZ File
	    TinkerXYZOutputFile xyzOutput = new TinkerXYZOutputFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + "_minimized.xyz");
        
        // retrieve file that contains energy
        TinkerMinimizationLogFile output = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
        
        // delete files after normal termination
        try
            {
                File[] files = new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles();
                for ( File f : files )
                    {
                        String filename = f.getName();
                        if ( f.getName().startsWith(baseFilename) )
                            {
                                f.delete();
                                //System.out.println(f.getName() + " deleted.");
                            }
                    }
            }
        catch (Exception e)
            {
                System.out.println("Error while trying to delete files:");
                e.printStackTrace();
            }

        // construct and return result
        return new TinkerMinimizationResult(output, elapsedTime, xyzOutput);
    }

    /** A class representing output files from a minimization call. Contains the molecule and log files produced from minimize call */
    public static class TinkerMinimizationResult implements Result, Serializable
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

	    /** The Tinker Minimiization Log File that is produced by minimize containing the energy of the optimized molecule */
	    public final TinkerMinimizationLogFile tinkerMinimizationLogFile;

	    /** The run time of the call in seconds. */
	    public final double elapsedTime;

	    /** The xyz file that is produced by a minimize call */
	    public final TinkerXYZOutputFile tinkerXYZOutputFile;

	    public TinkerMinimizationResult(TinkerMinimizationLogFile tinkerMinimizationLogFile, double elapsedTime, TinkerXYZOutputFile tinkerXYZOutputFile)
        {
	        this.tinkerMinimizationLogFile = tinkerMinimizationLogFile;
	        this.elapsedTime = elapsedTime;
	        this.tinkerXYZOutputFile = tinkerXYZOutputFile;
	    }
	
	    public Molecule getMolecule()
        {
	        return tinkerXYZOutputFile.molecule;
	    }

	    public double getEnergy()
        {
	        return tinkerMinimizationLogFile.energy;
	    }
	
	    @Override
	    public int hashCode()
        {
            return Objects.hash(tinkerMinimizationLogFile, elapsedTime, tinkerXYZOutputFile);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof TinkerMinimizationResult) )
                return false;

            TinkerMinimizationResult result = (TinkerMinimizationResult)obj;
            if ( Objects.equals(this.tinkerMinimizationLogFile, result.tinkerMinimizationLogFile) &&
                                this.elapsedTime == result.elapsedTime &&
                 Objects.equals(tinkerXYZOutputFile, result.tinkerXYZOutputFile) )
                return true;
            return false;
        }

        @Override
        public String toString()
        {
            return tinkerMinimizationLogFile.toString() + String.format(" time = %.3f s\n\n", elapsedTime) + tinkerXYZOutputFile.toString();
        }
    }
    
    /** A class representing the output from a minimization TINKER call */
    public static class TinkerMinimizationLogFile extends OutputFileFormat implements FileFormat
    {
	    /** The minimized energy in kcal/mol. */
	    public final Double energy;

        /** The final RMS gradient */
        public final Double gradient;

        /** The number of iterations */
        public final Integer iterations;

        /** the filename */
        public final String filename;

	    public TinkerMinimizationLogFile(String filename)
        {
            super(filename);
            this.filename = filename;

            // determine minimized energy
            Double readEnergy = null;
            Double readGradient = null;
            Integer readIterations = null;
            /*for ( int i=0; i < stringRepresentation.length(); i++ )
                {
                    String currentLetter = stringRepresentation.substring(i,i+1);
                    char thisChar = currentLetter.charAt(0);
                    int code = (int)thisChar;
                    System.out.println(currentLetter + " (" + code + ")");
                }*/
            for (String line : stringRepresentation.split("\n"))
                {
                    //System.out.println(index + " " + line);
                    String[] fields = line.trim().split("\\s+");
                    if ( line.indexOf("Final Function Value") > -1 )
                        readEnergy = Double.valueOf(fields[4]);
                    else if ( line.indexOf("Final RMS Gradient") > -1 )
                        {
                            readGradient = Double.valueOf(fields[4]);
                            //System.out.println(readGradient + " : " + line);
                        }
                    else 
                        {
                            try
                                {
                                    readIterations = Integer.valueOf(fields[0]);
                                }
                            catch (NumberFormatException e)
                                {
                                }
                        }
                }

            if ( readEnergy == null )
                throw new IllegalArgumentException("Energy not found!");
            if ( readGradient == null )
                throw new IllegalArgumentException("Gradient not found!");
            if ( readIterations == null )
                throw new IllegalArgumentException("Number of iterations not found!");
            this.energy = readEnergy;
            this.gradient = readGradient;
            this.iterations = readIterations;
        }

	    @Override
	    public int hashCode()
        {
            return Objects.hash(stringRepresentation, energy);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof TinkerMinimizationLogFile) )
                return false;

            TinkerMinimizationLogFile file = (TinkerMinimizationLogFile)obj;
            if ( Objects.equals(this.stringRepresentation, file.stringRepresentation) &&
                 Objects.equals(this.energy, file.energy ) )
                return true;
            return false;
        }

        @Override
        public String toString()
        {
            return String.format("TinkerMinimizationLogFile: energy = %.3f", energy);
        }
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(tinkerXYZInputFile, tinkerKeyFile);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof TinkerMinimizationJob) )
            return false;

        TinkerMinimizationJob job = (TinkerMinimizationJob)obj;
        if ( Objects.equals(this.tinkerXYZInputFile, job.tinkerXYZInputFile) &&
             Objects.equals(this.tinkerKeyFile,      job.tinkerKeyFile) )
            return true;
        return false;
    }

    @Override
    public String toString()
    {
        return "TinkerJob:\n" + tinkerXYZInputFile.stringRepresentation + "\n" + tinkerKeyFile.stringRepresentation;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        TinkerMinimizationLogFile file = new TinkerMinimizationLogFile("test.out");
        System.out.println(file.energy);
    }
}
