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
        List<Catalyst> conformations = new ArrayList<>();
        List<Double> energies = new ArrayList<>();
        List<List<Double>> distances = new ArrayList<>();
        List<String> distanceIdentifiers = new ArrayList<>();

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

        // read conformations
        int counter = 0;
        for ( List<Vector3D> geom : m.geometries )
            {
                // map the new conformation to the catalyst and print out gjfs
                Map<Atom,Atom> atomMap = new HashMap<>();
                if ( c.contents.size() != geom.size())
                    throw new IllegalArgumentException("Sizes of molecules do not match!");
                for ( int i = 0; i < c.contents.size(); i++ )
                    atomMap.put(c.getAtom(i+1), c.getAtom(i+1).moveAtom(geom.get(i)));
                c = c.moveAtoms(atomMap).setName(String.format("%05d", ++counter));

                conformations.add(c);
                String moleculeName = filename + String.format("%05d", counter);
                GaussianInputFile file = new GaussianInputFile(c, c.name, 
                    "#p opt b3lyp m06-2x geom=connect scrf=(solvent=benzene) freq");
                file.write(Settings.WORKING_DIRECTORY + "output/" + filename + "/" + moleculeName + ".gjf");
            }

        // read energies
        energies = new ArrayList<Double>(m.energies);

        // read inter(thio)urea distances
        // first, find pairs of (thio)ureas
        // then, we add an identifier to distanceIdentifiers
        // and all the distances are added to the list
        List<Fragment> fragments = c.fragmentList;
        for ( int i = 0; i < fragments.size(); i++)
            for ( int j = i; j < fragments.size(); j++)
                if ( !( fragments.get(i).ureaCarbon.element.equals(Element.DUMMY)
                        || fragments.get(j).ureaCarbon.element.equals(Element.DUMMY)) && (i!=j) )
                    {
                        distanceIdentifiers.add("C"+ c.getAtomNumber(fragments.get(i).ureaCarbon) + "-C" 
                                            + c.getAtomNumber(fragments.get(j).ureaCarbon) + " distance");
                        
                        List<Double> tempDistances = new ArrayList<>();
                        for ( Catalyst cat : conformations )
                            tempDistances.add(Vector3D.distance(
                                        cat.fragmentList.get(i).ureaCarbon.position,
                                        cat.fragmentList.get(j).ureaCarbon.position));
                    
                        distances.add(tempDistances);
                    }

        // copies of unsorted lists for reference
        final List<Catalyst> conformations2 = conformations;
        final List<List<Double>> distances2 = distances;
        final List<Double> energies2 = energies;
        
        // sort by energies
        Collections.sort(conformations,
            new Comparator<Catalyst>(){
                public int compare(Catalyst c1, Catalyst c2){
                    return Double.compare(energies2.get(conformations2.indexOf(c1)),
                                           energies2.get(conformations2.indexOf(c2)));
                }
            });

        for ( int n = 0; n < distances.size(); n++ )
        {
            final int n2 = n; // you didn't see this
            Collections.sort( distances.get(n), 
                new Comparator<Double>(){
                    public int compare(Double d1, Double d2){
                        return Double.compare(energies2.get(distances2.get(n2).indexOf(d1)),
                                                energies2.get(distances2.get(n2).indexOf(d2)));
                    }
                });
        }

        Collections.sort(energies);

        // write to file
        String fileString = "Conformation Number\tEnergy\t";
        for ( String s : distanceIdentifiers)
            fileString = fileString + "\t" + s;

        for ( int i = 0; i < energies.size(); i++ )
        {
            fileString = fileString + "\n";
            fileString = fileString + conformations.get(i).name + "\t\t\t\t";
            fileString = fileString + String.format("%10.6f", energies.get(i)) + "\t";
            for ( int j = 0; j < distances.size(); j++ )
                fileString = fileString + String.format("%-10.6f", distances.get(j).get(i)) + "\t";
        }
        System.out.println(fileString);
        try {
            PrintWriter writer = new PrintWriter(Settings.WORKING_DIRECTORY + "/output/" + filename + "/energy");
            writer.println(fileString);
            writer.close();
        }
        catch(IOException e)
        {
            e.printStackTrace();
        }

        // sort by distances
        for ( int n = 0; n < distances.size(); n++)
        {
            final int n2 = n; // somewhat suspicious, but please ignore
            Collections.sort(conformations,
                new Comparator<Catalyst>(){
                    public int compare(Catalyst c1, Catalyst c2){
                        return Double.compare(distances2.get(n2).get(conformations2.indexOf(c1)),
                                              distances2.get(n2).get(conformations2.indexOf(c2)));
                    }
                });
    
            Collections.sort(energies, 
                new Comparator<Double>(){
                    public int compare(Double d1, Double d2){
                        return Double.compare(distances2.get(n2).get(energies2.indexOf(d1)),
                                                    distances2.get(n2).get(energies2.indexOf(d2)));
                    }
                });

            for ( int m1 = 0; m1 < distances.size(); m1++)
            {
                final int m2 = m1;
                Collections.sort(distances.get(m1),
                    new Comparator<Double>(){
                        public int compare(Double d1, Double d2){
                            return Double.compare(distances2.get(n2).get(distances2.get(m2).indexOf(d1)),
                                                    distances2.get(n2).get(distances2.get(m2).indexOf(d2)));
                        }
                    });
            }

            // write to file
            fileString = "Conformation Number\tEnergy\t";
            for ( String s : distanceIdentifiers)
                fileString = fileString + "\t" + s;

            for ( int i = 0; i < energies.size(); i++ )
            {
                fileString = fileString + "\n";
                fileString = fileString + conformations.get(i).name + "\t\t\t\t";
                fileString = fileString + String.format("%10.6f", energies.get(i)) + "\t";
                for ( int j = 0; j < distances.size(); j++ )
                fileString = fileString + String.format("%-10.6f", distances.get(j).get(i)) + "\t";
            }
            System.out.println(fileString);
            try {
                PrintWriter writer = new PrintWriter(new FileOutputStream(new File(Settings.WORKING_DIRECTORY + "/output/" + filename + "/analysis"), true));
                writer.println(fileString);
                writer.close();
            }
            catch(IOException e)
            {
                e.printStackTrace();
            }
        }
    }              
}
