import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import Jama.*;

/**
* Represents a molecular fragment.  This class is effectively immutable and serializable.
* The connectivity graph is not exposed because it is potentially mutable.
*/
public class Fragment extends Molecule implements Immutable, Serializable
{
    public static final long serialVersionUID = 1L;

    /** The left connecting Atom of the fragment. */
    public final Atom leftConnect;

    /** The right connecting Atom of the fragment. */
    public final Atom rightConnect;

    /** The Atom that is the (thio)urea carbon, if any in this fragment. */
    public final Atom ureaCarbon;

    /** The FragmentType of this fragment. */
    public FragmentType fragmentType;

    /** 
    * Factory method to construct Fragment from GJFFragment. 
    * @param fragmentFile a GJFfragment read from a .gjf metafile
    * @return a Fragment
    */
    public static Fragment createFragment(GJFfragment fragmentFile)
    {
        Molecule molecule = fragmentFile.molecule;
        String name = molecule.name;
        List<Atom> contents = molecule.contents;
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = molecule.connectivity;
        Atom leftConnect = fragmentFile.leftConnect;
        Atom rightConnect = fragmentFile.rightConnect;
        Atom ureaCarbon = fragmentFile.ureaCarbon;
        FragmentType fragmentType = fragmentFile.fragmentType;

        return new Fragment(name, contents, connectivity, leftConnect, rightConnect, ureaCarbon, fragmentType);
    }

    /**
    * Private constructor for Fragment.
    * @param name the name of the Fragment
    * @param contents a list of Atoms in the Fragment
    * @param connectivity a graph of the connectivity
    * @param leftConnect
    * @param rightConnect
    * @param ureaCarbon the (thio)urea carbon, if any
    * @param fragmentType
    * @return a Fragment
    */

    private Fragment(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, Atom leftConnect, Atom rightConnect, Atom ureaCarbon, FragmentType fragmentType)
    {
        super(name, contents, connectivity);
        this.leftConnect = leftConnect;
        this.rightConnect = rightConnect;
        this.ureaCarbon = ureaCarbon;
        this.fragmentType = fragmentType;
    }

    /**
    * Returns a String representation of the Fragment.
    * @return a String
    */
    @Override
    public String toString()
    {
        String superToString = super.toString();
        superToString += "\nLeft: " + leftConnect + "\nRight " + rightConnect;
        superToString += "\nUrea Carbon: " + ureaCarbon + "\nfragmentType: " + fragmentType;
        return superToString;
    }

    /**
    * Considers whether a given Object is equal to this one.
    * @param obj the object
    * @return true if equal
    */
    @Override
    public boolean equals(Object obj)
    {
        boolean superEquals = super.equals(obj);
        Fragment fragment = (Fragment)obj;
        return superEquals && (fragment.leftConnect == this.leftConnect) && (fragment.rightConnect == this.rightConnect) && (fragment.fragmentType == this.fragmentType) && (fragment.ureaCarbon == this.ureaCarbon);
    }

    /**
    * Returns the hash code of this Fragment/
    * @return the hash code
    */
    @Override
    public int hashCode()
    {
        return Objects.hash(name, contents, connectivity, leftConnect, rightConnect, ureaCarbon, fragmentType);
    }

}

