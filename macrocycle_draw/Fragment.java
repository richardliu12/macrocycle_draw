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
     * Factory method to create a Fragment given a map of old atoms to new atoms. 
     * Should be used to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     * @return a Fragment
     */
    public Fragment moveAtoms(Map<Atom,Atom> atomMap)
    {
        Molecule molecule = new Molecule(name, contents, connectivity);
        Molecule newMolecule = molecule.moveAtoms(atomMap);
        Atom newLeftConnect = this.leftConnect;
        Atom newRightConnect = this.rightConnect;
        Atom newUreaCarbon = this.ureaCarbon;
        
        List<Atom> newContents = new LinkedList<Atom>();
        if ( atomMap.containsKey(this.leftConnect) )
                    newLeftConnect = atomMap.get(this.leftConnect);
        if ( atomMap.containsKey(this.rightConnect) )
                    newRightConnect = atomMap.get(this.rightConnect);
        if ( atomMap.containsKey(this.ureaCarbon) )
                    newUreaCarbon = atomMap.get(this.ureaCarbon);

        return new Fragment(newMolecule.name, newMolecule.contents, newMolecule.connectivity, newLeftConnect, newRightConnect, newUreaCarbon, this.fragmentType);
    }

    /**
     * Factory method to create a new Fragment by transforming this one.
     * The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation that we apply to the Atoms in this to get the Atoms in the output
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this rotated by rot, plus shift
     */
    public Fragment transform(Rotation rot, Vector3D shift)
    {
        Molecule molecule = new Molecule(name, contents, connectivity);
        Molecule newMolecule = molecule.transform(rot, shift);
        Atom newLeftConnect = this.leftConnect.transform(rot,shift);
        Atom newRightConnect = this.rightConnect.transform(rot,shift);
        Atom newUreaCarbon = this.ureaCarbon.transform(rot,shift);
        
        return new Fragment(newMolecule.name, newMolecule.contents, newMolecule.connectivity, newLeftConnect, newRightConnect, newUreaCarbon, this.fragmentType);
    }

    /**
     * Factory method to create a new Fragment by shifting this one.
     * Uses Rotation.IDENTITY as a parameter in the transform method
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this plus shift
     */
    public Fragment shift(Vector3D shift)
    {
        return transform(Rotation.IDENTITY, shift);
    }

    /**
     * Factory method to create a new Fragment by shifting this one to have its center at the origin.
     * @return this shifted by minus the barycenter of this 
     */
    public Fragment normalize()
    {
        Vector3D barycenter = this.getCentroid();
        return shift(barycenter.negate());
    }

    /**
    * Factory method to create a new Fragment by modifying the type of this one
    */
    public Fragment newFragmentType(FragmentType ftype)
    {
        return new Fragment(name, contents, connectivity, leftConnect, rightConnect, ureaCarbon, ftype);
    }

    /**
     * Moves the group associated with atom2 to the specified distance.
     * Motion occurs along the atom1-atom2 bond vector.  Note that this returns a new
     * molecule.  No checks are made.
     * @param atom1 this atom will be held fixed
     * @param atom2 this atom and anything connected to it will be moved
     * @param requestedDistance the requested distance in Angstroms
     * @return a new Fragment containing the same connectivity but new positions
     */
    public Fragment setDistance(Atom atom1, Atom atom2, double requestedDistance)
    {
       Map<Atom,Atom> atomMap = setDistanceMap(atom1, atom2, requestedDistance);
       return moveAtoms(atomMap); 
    }

    /**
     * Alias method.  Atom indices are 1, 2, ..., n.  No checks.
     */
    public Fragment setDistance(int i, int j, double requestedDistance)
    {
        return setDistance(contents.get(i-1), contents.get(j-1), requestedDistance);
    }

    /**
     * Rotates the atom1-atom2-atom3 angle, moving only atom3 and anything in its
     * attached subgraph.  No checks.
     * @param atom1 will not be moved
     * @param atom2 will not be moved
     * @param atom3 will be moved
     * @param theta rotation in degrees
     * @return this rotated
     */
    public Fragment rotateAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        Map<Atom,Atom> atomMap2 = rotateAngleMap(atom1, atom2, atom3, theta);
        return moveAtoms(atomMap2);
    }
   
    /**
     * Method alias.  Indices are 1,2,...,n.  No checks.
     */
    public Fragment rotateAngle(int i, int j, int k, double theta)
    {
        return rotateAngle(contents.get(i-1), contents.get(j-1), contents.get(k-1), theta);
    }

    /**
     * Set the atom1-atom2-atom3 angle to theta degrees, moving atom3 and its subgraph only.
     * New Fragment returned.  No checks.
     * @param atom1 not moved
     * @param atom2 not moved
     * @param atom3 moved
     * @param theta desired angle in degrees
     */
    public Fragment setAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        double currentAngle = getAngle(atom1, atom2, atom3);
        double requiredRotation = theta - currentAngle;
        return rotateAngle(atom1, atom2, atom3, requiredRotation);
    }

    /**
     * Method alias. 
     */
    public Fragment setAngle(int i, int j, int k, double theta)
    {
        return setAngle(contents.get(i-1), contents.get(j-1), contents.get(k-1), theta);
    }

    /**
     * Returns a new Fragment with a rotated dihedral.
     * Note that the old AtomTorsion will no longer point to the new Molecule.
     * @param theta the desired dihedral angle in degrees
     * @return the new Fragment
     */
    public Fragment setDihedral(AtomTorsion atomTorsion, double theta)
    {
        Map<Atom,Atom> atomMap2 = setDihedralMap(atomTorsion, theta);
        return moveAtoms(atomMap2);
    }

    public Fragment setDihedral(ProtoTorsion protoTorsion, double theta)
    {
        AtomTorsion atomTorsion = protoTorsion.getAtomTorsion(this);
        return setDihedral(atomTorsion, theta);
    }

    /**
     * Returns a new Molecule with a rotated dihedral.  Alias method.
     * Note that the old IndexTorsion will still be valid for the new Molecule.
     */
    public Molecule setDihedral(IndexTorsion indexTorsion, double theta)
    {
        return setDihedral(indexTorsion.getAtomTorsion(this), theta);
    }

    /**
     * Creates a new Fragment where atom2 and its subgraph have been moved to
     * make atom1 sp2-hybridized (bond angles set at 120 degrees).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @param forceAngle true if we want to force the atom1alpha-atom1-atom1beta angle to 120 (safe for non-prolines)
     * @return a new Fragment with adjusted hybridization
     */
    public Fragment set_sp2(Atom atom1, Atom atom2, boolean forceAngle)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        Map<Atom,Atom> newAtomMap = set_sp2_map(atom1, atom2, forceAngle);
        Fragment rotatedFragment = moveAtoms(newAtomMap);

        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Fragment returnFragment = rotatedFragment.setDistance(atom1number, atom2number, currentLength); 
        return returnFragment;
    }

    /**
    * Alias method.  
    */
    public Fragment set_sp2(Atom atom1, Atom atom2)
    {
        return set_sp2(atom1, atom2, true);
    }

    /**
     * Creates a new Fragment where atom2 and its subgraph have been moved to
     * make atom1 sp3-hybridized (tetrahedral geometry).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @return a new Fragment with adjusted hybridization
     */
    public Fragment set_sp3(Atom atom1, Atom atom2)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        Map<Atom,Atom> newAtomMap = set_sp3_map(atom1, atom2);
        Fragment rotatedFragment = moveAtoms(newAtomMap);

        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Fragment returnFragment = rotatedFragment.setDistance(atom1number, atom2number, currentLength); 
        return returnFragment;
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
      	if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Fragment) )
            return false;

        Molecule anotherMolecule = (Molecule)obj;
        if ( this.name.equals(anotherMolecule.name) &&
             this.contents.equals(anotherMolecule.contents) &&
             this.connectivity.equals(anotherMolecule.connectivity) )
        {   
		Fragment fragment = (Fragment)obj;
        	return Objects.equals(fragment.leftConnect, this.leftConnect) && Objects.equals(fragment.rightConnect, this.rightConnect) && Objects.equals(fragment.fragmentType, this.fragmentType) && Objects.equals(fragment.ureaCarbon, this.ureaCarbon);
    	}
	
	return false;
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

