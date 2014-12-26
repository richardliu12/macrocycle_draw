import java.io.*;
import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/**
 * Represents a Gaussian input file with metadata to describe a fragment.
  */
public class GJFfragment extends OutputFileFormat implements Immutable
{
    /** The Molecule representation of this file. */
    public final Molecule molecule;
    
    /** The left connecting Atom of the fragment. */
    public final Atom leftConnect;

    /** The right connecting Atom of the fragment. */
    public final Atom rightConnect;

    /** The Atom that is the (thio)urea carbon, if any in this fragment. */
    public final Atom ureaCarbon;

    /** The FragmentType of this fragment. */
    public FragmentType fragmentType;

    /** A list of chiral atoms. */
    public final ImmutableList<Atom> chiralAtoms;

    /** A graph of rotatable bonds. */
    protected SimpleWeightedGraph<Atom,DefaultWeightedEdge> rotatableBonds;
 
    /**
     * Reads the geometry and connectivity, as well as metadata.
     * @param filename the location of the gjf file
    */
    public GJFfragment(String filename)
    {
	    super(filename);
	
        // read geometry
	    String name = filename.split("\\.")[0];
        List<Atom> contents = new ArrayList<>();
	    SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        int blanks = 0;
        boolean lastBlank = false;
        boolean inGeometryBlock = false;
        int leftConnectAtomNumber = 0;
        int rightConnectAtomNumber = 0;
        int ureaCarbonAtomNumber = 0;
        List<Integer> chiralAtomNumbers = new ArrayList<>();
        SimpleWeightedGraph<Integer,DefaultWeightedEdge> rotatableBondIndices = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);

        for (List<String> line : fileContents)
            {
                // keep track of how many blanks we have seen
                if ( line.size() == 1 && line.get(0).length() == 0 )
                    {
                        if ( lastBlank == false )
                            {
                                blanks++;
                                lastBlank = true;
                            }
                        continue;
                    }
                else
                    lastBlank = false;

                // read the metadata 
                if ( blanks == 1 )
                    {   
                        for (String s : line)
                        {  
                            if ( s.split("@").length <= 2 )
                                continue;
                                
                            if ( s.split("@")[1].toLowerCase().equals("left_connect") )
                            {
                                if ( s.split("@").length != 3 )
                                    throw new IllegalArgumentException("improper atom number specification for left_connect in " + filename + ":\n" + line.toString());
                                else
                                    leftConnectAtomNumber = Integer.parseInt(s.split("@")[2]);
                            }
                            else if ( s.split("@")[1].toLowerCase().equals("right_connect") )
                            {
                                if ( s.split("@").length != 3 )
                                    throw new IllegalArgumentException("improper atom number specification for right_connect in " + filename + ":\n" + line.toString());
                                else
                                    rightConnectAtomNumber = Integer.parseInt(s.split("@")[2]);
                            }
                            else if ( s.split("@")[1].toLowerCase().equals("urea_carbon") )
                            {
                                if ( s.split("@").length != 3 )
                                    throw new IllegalArgumentException("improper atom number specification for urea_carbon in " + filename + ":\n" + line.toString());
                                else
                                    ureaCarbonAtomNumber = Integer.parseInt(s.split("@")[2]);
                            }
                            else if ( s.split("@")[1].toLowerCase().equals("fragment_type") )
                            {
                                if ( s.split("@").length != 3 )
                                    throw new IllegalArgumentException("improper fragment_type in " + filename + ":\n" + line.toString());
                                else switch( s.split("@")[2] ){
                                    case "urea": fragmentType = FragmentType.UREA;
                                        break;
                                    case "thiourea": fragmentType = FragmentType.THIOUREA;
                                        break;
                                    case "linker_1": fragmentType = FragmentType.LINKER_1;
                                        break;
                                    case "linker_2": fragmentType = FragmentType.LINKER_2;
                                        break;
                                    case "linker_3": fragmentType = FragmentType.LINKER_3;
                                        break;
                                    case "linker_4": fragmentType = FragmentType.LINKER_4;
                                        break;
                                    default: throw new IllegalArgumentException("improper fragment_type in " + filename + ":\n" + line.toString());
                                }
                            }
                            else if ( s.split("@")[1].toLowerCase().equals("chiral_atom") )
                            {
                               if ( s.split("@").length < 3 ) 
                                   throw new IllegalArgumentException("improper chiral_atom specification in " + filename + ":\n" + line.toString());
                                else
                                    for ( int i = 2; i < s.split("@").length; i++)
                                        chiralAtomNumbers.add(Integer.parseInt(s.split("@")[i]));
                            }
                            else if ( s.split("@")[1].toLowerCase().equals("rotatable_bond") )
                            {
                               if ( s.split("@").length < 4 )
                                   throw new IllegalArgumentException("improper rotatable_bond specification in " + filename + ":\n" + line.toString());
                                else 
                                    for ( int i = 2; i+1 < s.split("@").length; i++)
                                    {
                                        int atomNumber1 = Integer.parseInt(s.split("@")[i])-1;
                                        int atomNumber2 = Integer.parseInt(s.split("@")[i+1])-1;
                                        
                                        rotatableBondIndices.addVertex(atomNumber1);
                                        rotatableBondIndices.addVertex(atomNumber2);
                                        rotatableBondIndices.addEdge(atomNumber1,atomNumber2);
                                    }
                            }
                        }
                        continue;
                    }
                else if ( blanks != 2 )
                    continue;

                // deal with the charge and multiplicity card (by ignoring it)
                if ( line.size() == 2 && inGeometryBlock == false )
                    {
                        inGeometryBlock = true;
                        continue;
                    }

                if ( line.size() != 4 && inGeometryBlock == false )
                    throw new IllegalArgumentException("unexpected text in geometry block in " + filename + ":\n" + line.toString());
                
                // create atom
                // tinker atom types will be nonsense, of course
                Atom newAtom = new Atom(line.get(0), new Vector3D(Double.parseDouble(line.get(1)), Double.parseDouble(line.get(2)), Double.parseDouble(line.get(3))), 1);
                contents.add(newAtom);
                connectivity.addVertex(newAtom);
            }

        // read connectivity
        blanks = 0;
        lastBlank = false;
        for (List<String> line : fileContents)
            {
                // read the fourth block of text
                if ( line.size() == 1 && line.get(0).length() == 0 )
                    {
                        if ( lastBlank == false )
                            {
                                blanks++;
                                lastBlank = true;
                            }
                        continue;
                    }
                else
                    lastBlank = false;
               
               // only read connectivity lines
                if ( blanks != 3 )
                    continue;

                Atom fromAtom = contents.get(Integer.parseInt(line.get(0))-1);
                for (int i=1; i < line.size(); i+=2)
                    {
                        int toAtomIndex = Integer.parseInt(line.get(i))-1;
                        Atom toAtom = contents.get(toAtomIndex);
                        double bondOrder = Double.parseDouble(line.get(i+1));
                        DefaultWeightedEdge thisEdge = connectivity.addEdge(fromAtom, toAtom);
                        connectivity.setEdgeWeight(thisEdge, bondOrder);
                    }
            }

        // create the molecule, label atoms with metadata
	    molecule = new Molecule(name, contents, connectivity);
        leftConnect = molecule.getAtom(leftConnectAtomNumber);
        rightConnect = molecule.getAtom(rightConnectAtomNumber);
        if (ureaCarbonAtomNumber != 0)
            ureaCarbon = molecule.getAtom(ureaCarbonAtomNumber);
        else
            ureaCarbon = new Atom("Q", new Vector3D(0,0,0), 0);
        
        List<Atom> tempChiralAtoms = new ArrayList<>();
        for ( int i : chiralAtomNumbers )
            tempChiralAtoms.add(molecule.getAtom(i));
        chiralAtoms = ImmutableList.copyOf(tempChiralAtoms);
        
        rotatableBonds = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        for ( int atomNumber : rotatableBondIndices.vertexSet() )
            rotatableBonds.addVertex(contents.get(atomNumber));
        for ( DefaultWeightedEdge e : rotatableBondIndices.edgeSet() )
            {
                int atomNumber1 = rotatableBondIndices.getEdgeSource(e);
                int atomNumber2 = rotatableBondIndices.getEdgeTarget(e);
                rotatableBonds.addEdge(contents.get(atomNumber1), contents.get(atomNumber2));
            }                    
    }
        /**
        * For testing.
        */
        public static void main(String[] args)
        {
            GJFfragment g1 = new GJFfragment("test.gjf");
            Fragment f1 = Fragment.createFragment(g1);
            GJFfragment g2 = new GJFfragment("test2.gjf");
            Fragment f2 = Fragment.createFragment(g2);
            GJFfragment g3 = new GJFfragment("test3.gjf");
            Fragment f3 = Fragment.createFragment(g3);

            Catalyst cat = new Catalyst(f1);
            cat = cat.addLeft(f2);
            cat = cat.addRight(f3);
            System.out.println(cat.getTorsions());
            System.out.println(cat.getOPLSenergy());

            Molecule cyc = MonteCarlo.cyclize(cat, cat.getTorsions(), cat.getAtomNumber(cat.getLeftConnect()), cat.getAtomNumber(cat.getRightConnect()));
            MOL2InputFile m = new MOL2InputFile(cyc);
            m.write("test_join.mol2");

        }
}
