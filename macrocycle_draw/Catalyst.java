import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import Jama.*;

/**
* Represents a collection of connected fragments.  This class is effectively immutable
* and serializable.
*/
public class Catalyst extends Molecule implements Immutable, Serializable
{
    public static final long serialVersionUID = 1L;

    // the default bond length for newly formed bonds between fragments
    final double BOND_LENGTH = 1.35;

    /** Ordered list of all Fragments in this Catalyst */
    public final ImmutableList<Fragment> fragmentList;

    /**
    * Private constructor used to create modified copies of the Catalyst.
    * External modification to the Catalyst should be performed using methods
    * such as addLeft() and addRight() only.
    * @param name the name of the catalyst
    * @param contents an List of the Atoms
    * @param connectivity a graph of the bonds
    * @param fragmentList a List of the Fragments contained
    */
    protected Catalyst(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, List<Fragment> fragmentList)
    {
        super(name, contents, connectivity);
        this.fragmentList = ImmutableList.copyOf(fragmentList);
    }

    /**
    * Returns a deep copy with new name.
    */
    @Override
    public Catalyst setName(String name)
    {
        return new Catalyst(name, this.contents, this.connectivity, this.fragmentList);
    }

    /**
    * Default constructor that initializes the Catalyst with the first fragment 
    * @param fragment the first fragment
    */
    public Catalyst(Fragment fragment)
    {
        super(fragment.name, fragment.contents, fragment.connectivity);
        List<Fragment> newFragmentList = new ArrayList<Fragment>();
        newFragmentList.add(fragment);
        this.fragmentList = ImmutableList.copyOf(newFragmentList);
    }
   
    /**
    * Returns the left connection point of the linear chain.
    */
    public Atom getLeftConnect()
    {
        return fragmentList.get(0).leftConnect;
    }

    public Atom getRightConnect()
    {
        return fragmentList.get(fragmentList.size() - 1).rightConnect;
    }

    /**
     * Generates a list of IndexTorsions from the rotatable bonds.
     * @return a list of IndexTorsions
     */
    public List<IndexTorsion> getTorsions()
    {
        List<IndexTorsion> returnTorsions = new ArrayList<>();

        for ( Fragment f : fragmentList )
        {
            for ( DefaultWeightedEdge e : f.rotatableBonds.edgeSet() )
                {
                    Atom fromAtom = f.rotatableBonds.getEdgeSource(e);
                    Atom toAtom = f.rotatableBonds.getEdgeTarget(e);

                    List<Atom> fromAtomNeighbors = new ArrayList<Atom>(this.getAdjacentAtoms(fromAtom));
                    List<Atom> toAtomNeighbors = new ArrayList<Atom>(this.getAdjacentAtoms(toAtom));

                    fromAtomNeighbors.remove(toAtom);
                    toAtomNeighbors.remove(fromAtom);

                    if ( fromAtomNeighbors.size() == 0 || toAtomNeighbors.size() == 0 )
                        continue;
                    else
                        returnTorsions.add(IndexTorsion.createIndexTorsion(getAtomNumber(fromAtomNeighbors.get(0)), getAtomNumber(fromAtom), getAtomNumber(toAtom), getAtomNumber(toAtomNeighbors.get(0)), this));
                }
            if ( fragmentList.indexOf(f) + 1 < fragmentList.size() )
                {
                    Atom fromAtom = f.rightConnect;
                    Atom toAtom = fragmentList.get(fragmentList.indexOf(f) + 1).leftConnect;
                    
                    List<Atom> fromAtomNeighbors = new ArrayList<Atom>(this.getAdjacentAtoms(fromAtom));
                    List<Atom> toAtomNeighbors = new ArrayList<Atom>(this.getAdjacentAtoms(toAtom));

                    fromAtomNeighbors.remove(toAtom);
                    toAtomNeighbors.remove(fromAtom);

                    if ( fromAtomNeighbors.size() != 0 && toAtomNeighbors.size() != 0 )
                        returnTorsions.add(IndexTorsion.createIndexTorsion(getAtomNumber(fromAtomNeighbors.get(0)), getAtomNumber(fromAtom), getAtomNumber(toAtom), getAtomNumber(toAtomNeighbors.get(0)), this));
                }
        }
        return returnTorsions;
    }
    
    /**
     * Gets torsions from a cyclic catalyst.  Cuts the connection bond and returns the
     * remaining torsions.
     */
     public List<IndexTorsion> getLinearTorsions()
     {
         this.connectivity.removeEdge(fragmentList.get(0).leftConnect, fragmentList.get(fragmentList.size()-1).rightConnect);
         List<IndexTorsion> returnTorsions = getTorsions();
         this.connectivity.addEdge(fragmentList.get(0).leftConnect, fragmentList.get(fragmentList.size()-1).rightConnect);
         return returnTorsions;
     }

    /**
     * Cyclization method.  Returns a cyclized version of this catalyst.
     */
     public Catalyst cyclize()
     {
         Molecule m = MonteCarlo.cyclize(this, this.getTorsions(), getAtomNumber(this.getLeftConnect()), getAtomNumber(this.getRightConnect()));
         Map<Atom,Atom> atomMap = this.matchMap(m);
         // create the new bond
         DefaultWeightedEdge e = connectivity.addEdge(this.getLeftConnect(), this.getRightConnect());
         connectivity.setEdgeWeight(e, 1);

         return this.moveAtoms(atomMap);
     }

    /** 
     * Factory method that returns a copy of the Catalyst with fragment appended on left.
     * The catalyst left connection atom is connected to the fragment's right connection
     * point, and angles are crudely adjusted.
     * @param fragment the fragment to be added
     * @return elongated Catalyst
     */
    @SuppressWarnings("unchecked")
    public Catalyst addLeft(Fragment fragment)
    {   
        if ( fragmentList.size() == 0 )
        {
            List<Fragment> newFragmentList = new ArrayList<Fragment>();
            newFragmentList.add(fragment);
            return new Catalyst(fragment.name, fragment.contents, fragment.connectivity, newFragmentList);
        }
        else
        {
            String newName = fragment.name + "_" + name;

            List<Atom> newContents = new ArrayList<Atom>(fragment.contents);
            newContents.addAll(contents);

            SimpleWeightedGraph<Atom, DefaultWeightedEdge> newC = (SimpleWeightedGraph<Atom, DefaultWeightedEdge>)connectivity.clone();
            Graphs.addGraph(newC,fragment.connectivity);
            Graphs.addEdge(newC, this.fragmentList.get(0).leftConnect,fragment.rightConnect, 1.0);

            List<Fragment> newFragmentList = new ArrayList<Fragment>();
            newFragmentList.add(fragment);
            newFragmentList.addAll(fragmentList);

            Catalyst returnCatalyst = new Catalyst(newName, newContents, newC, newFragmentList);
            return returnCatalyst.geometryCorrect(this.fragmentList.get(0).leftConnect, fragment.rightConnect);
        }
    }
    
    /** 
     * Factory method that returns a copy of the Catalyst with fragment appended on right.
     * The catalyst right connection atom is connected to the fragment's left connection
     * point, and angles are crudely adjusted.
     * @param fragment the fragment to be added
     * @return elongated Catalyst
     */
    @SuppressWarnings("unchecked") 
    public Catalyst addRight(Fragment fragment)
    {
        if ( fragmentList.size() == 0 )
        {
            List<Fragment> newFragmentList = new ArrayList<Fragment>();
            newFragmentList.add(fragment);
            return new Catalyst(fragment.name, fragment.contents, fragment.connectivity, newFragmentList);
        }
        else
        {
            String newName = name + "_" + fragment.name;

            List<Atom> newContents = new ArrayList<Atom>(contents);
            newContents.addAll(fragment.contents);

            SimpleWeightedGraph<Atom, DefaultWeightedEdge> newC = (SimpleWeightedGraph<Atom, DefaultWeightedEdge>)fragment.connectivity.clone();
            Graphs.addGraph(newC,connectivity);
            Graphs.addEdge(newC, this.fragmentList.get(this.fragmentList.size()-1).rightConnect,fragment.leftConnect, 1.0);

            List<Fragment> newFragmentList = new ArrayList<Fragment>();
            newFragmentList.addAll(fragmentList);
            newFragmentList.add(fragment);

            Catalyst returnCatalyst = new Catalyst(newName, newContents, newC, newFragmentList);
            return returnCatalyst.geometryCorrect(this.fragmentList.get(this.fragmentList.size()-1).rightConnect, fragment.leftConnect);
        }
    }

    /**
     * Corrects the geometry at the given atoms.
     * The geometry correction sets the bond length to a default value.
     * @param atom1 the first atom
     * @param atom2 the second atom
     */
     public Catalyst geometryCorrect(Atom atom1, Atom atom2)
     {
        int atomNumber1 = getAtomNumber(atom1);
        int atomNumber2 = getAtomNumber(atom2);

        // obtain the set of atoms adjacent to atom 1, not including atom 2
        Set<Atom> adjSet1 = getAdjacentAtoms(atom1);
        adjSet1.remove(atom2);
        List<Atom> adj1 = new LinkedList<Atom>(adjSet1);

    	Catalyst newCatalyst1;

        // determine how many bonds to atom 1, then correct geometry appropriately
        switch(adj1.size())
        {
            case 3: newCatalyst1 = set_sp3(atom1, atom2);
                break;
            case 2: newCatalyst1 = set_sp2(atom1, atom2);
                break;
            case 1: newCatalyst1 = setAngle(adj1.get(0), atom1, atom2, 120);
                break;
            default: throw new IllegalArgumentException("Atom has too many or too few bonds! (" + adj1.size() + ") \n" + atom1.toString());
        }

        // find moved atoms and repeat procedure for atom 2
        atom1 = newCatalyst1.getAtom(atomNumber1);
        atom2 = newCatalyst1.getAtom(atomNumber2);

        Set<Atom> adjSet2 = newCatalyst1.getAdjacentAtoms(atom2);
        adjSet2.remove(atom1);
        List<Atom> adj2 = new LinkedList<Atom>(adjSet2);

    	Catalyst newCatalyst2;

        switch(adj2.size())
        {
            case 3: newCatalyst2 = newCatalyst1.set_sp3(atom2, atom1);
                break;
            case 2: newCatalyst2 = newCatalyst1.set_sp2(atom2, atom1);
                break;
            case 1: newCatalyst2 = newCatalyst1.setAngle(adj2.get(0), atom2, atom1, 120);
                break;
   
            default: throw new IllegalArgumentException("Atom has too many or too few bonds! (" + adj2.size() +") \n" + atom1.toString());
        }
    
        // correct distance
        Catalyst returnCatalyst = newCatalyst2.setDistance(atomNumber1, atomNumber2, BOND_LENGTH);
        return returnCatalyst.normalize();
     }

    /**
     * Factory method to create a Catalyst given a map of old atoms to new atoms. 
     * Should be used to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     * @return a moved Catalyst
     */
    public Catalyst moveAtoms(Map<Atom,Atom> atomMap)
    {
        Molecule newMolecule = super.moveAtoms(atomMap);
        List<Fragment> newFragmentList = new ArrayList<Fragment>();
        for (Fragment fragment : fragmentList)
        {
            fragment = fragment.moveAtoms(atomMap);
            newFragmentList.add(fragment);
        }

        return new Catalyst(newMolecule.name, newMolecule.contents, newMolecule.connectivity, newFragmentList);
    }

    /**
     * Factory method to create a new Catalyst by transforming this one.
     * The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation that we apply to the Atoms in this to get the Atoms in the output
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this rotated by rot, plus shift
     */
    public Catalyst transform(Rotation rot, Vector3D shift)
    {
        Molecule newMolecule = super.transform(rot, shift);
        List<Fragment> newFragmentList = new ArrayList<Fragment>();
        for (Fragment fragment : fragmentList)
        {
            fragment = fragment.transform(rot,shift);
            newFragmentList.add(fragment);
        }

        return new Catalyst(newMolecule.name, newMolecule.contents, newMolecule.connectivity, newFragmentList);
    }

    /**
     * Factory method to create a new Catalyst by shifting this one.
     * Uses Rotation.IDENTITY as a parameter in the transform method
     * @param shift a vector that we add to the positions of the Atoms in this to get the positions of the Atoms in the output
     * @return this plus shift
     */
    public Catalyst shift(Vector3D shift)
    {
        return transform(Rotation.IDENTITY, shift);
    }

    /**
     * Moves the group associated with atom2 to the specified distance.
     * Motion occurs along the atom1-atom2 bond vector.  Note that this returns a new
     * molecule.  No checks are made.
     * @param atom1 this atom will be held fixed
     * @param atom2 this atom and anything connected to it will be moved
     * @param requestedDistance the requested distance in Angstroms
     * @return a new Catalyst containing the same connectivity but new positions
     */
    public Catalyst setDistance(Atom atom1, Atom atom2, double requestedDistance)
    {
       Map<Atom,Atom> atomMap = setDistanceMap(atom1, atom2, requestedDistance);
       return moveAtoms(atomMap); 
    }

    /**
     * Alias method.  Atom indices are 1, 2, ..., n.  No checks.
     */
    public Catalyst setDistance(int i, int j, double requestedDistance)
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
    public Catalyst rotateAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        Map<Atom,Atom> atomMap2 = rotateAngleMap(atom1, atom2, atom3, theta);
        return moveAtoms(atomMap2);
    }

    /**
     * Set the atom1-atom2-atom3 angle to theta degrees, moving atom3 and its subgraph only.
     * New Catalyst returned.  No checks.
     * @param atom1 not moved
     * @param atom2 not moved
     * @param atom3 moved
     * @param theta desired angle in degrees
     */
    public Catalyst setAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        double currentAngle = getAngle(atom1, atom2, atom3);
        double requiredRotation = theta - currentAngle;
        return rotateAngle(atom1, atom2, atom3, requiredRotation);
    }

    /**
     * Returns a new Catalyst with a rotated dihedral.
     * Note that the old AtomTorsion will no longer point to the new Molecule.
     * @param theta the desired dihedral angle in degrees
     * @return the new Catalyst
     */
    public Catalyst setDihedral(AtomTorsion atomTorsion, double theta)
    {
        Map<Atom,Atom> atomMap2 = setDihedralMap(atomTorsion, theta);
        return moveAtoms(atomMap2);
    }

    public Catalyst setDihedral(ProtoTorsion protoTorsion, double theta)
    {
        AtomTorsion atomTorsion = protoTorsion.getAtomTorsion(this);
        return setDihedral(atomTorsion, theta);
    }

    /**
     * Creates a new Catalyst where atom2 and its subgraph have been moved to
     * make atom1 sp2-hybridized (bond angles set at 120 degrees).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @param forceAngle true if we want to force the atom1alpha-atom1-atom1beta angle to 120 (safe for non-prolines)
     * @return a new Catalyst with adjusted hybridization
     */
    public Catalyst set_sp2(Atom atom1, Atom atom2, boolean forceAngle)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        Map<Atom,Atom> newAtomMap = set_sp2_map(atom1, atom2, forceAngle);
        Catalyst rotatedCatalyst = moveAtoms(newAtomMap);

        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        return rotatedCatalyst.setDistance(atom1number, atom2number, currentLength); 
    }

    /**
    * Alias method.  
    */
    public Catalyst set_sp2(Atom atom1, Atom atom2)
    {
        return set_sp2(atom1, atom2, true);
    }

    /**
     * Creates a new Catalyst where atom2 and its subgraph have been moved to
     * make atom1 sp3-hybridized (tetrahedral geometry).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @return a new Catalyst with adjusted hybridization
     */
    public Catalyst set_sp3(Atom atom1, Atom atom2)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);
        int atom1number = getAtomNumber(atom1);
        int atom2number = getAtomNumber(atom2);

        Map<Atom,Atom> newAtomMap = set_sp3_map(atom1, atom2);
        Catalyst rotatedCatalyst = moveAtoms(newAtomMap);

      	return rotatedCatalyst.setDistance(atom1number, atom2number, currentLength); 
    }

    /**
     * Factory method to create a new Catalyst by shifting this one to have its center at the origin.
     * @return this shifted by minus the barycenter of this 
     */
    public Catalyst normalize()
    {
        Vector3D barycenter = this.getCentroid();
        return shift(barycenter.negate());
    }

    /**
     * Returns the hash code of this Catalyst.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(name, contents, connectivity, fragmentList);
    }

    /**
     * Checks for object equality by comparing all fields.
     * @param obj the object to be compared to
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Catalyst) )
            return false;

        Catalyst anotherCatalyst = (Catalyst)obj;
        if ( this.name.equals(anotherCatalyst.name) &&
             this.contents.equals(anotherCatalyst.contents) &&
             this.connectivity.equals(anotherCatalyst.connectivity) && 
             this.fragmentList.equals(anotherCatalyst.fragmentList))
            return true;
        return false;
    }

    /**
     * For testing only.
     */
    public static void main(String[] args)
    {
        GJFfragment gjf = new GJFfragment("test.gjf");
        Fragment frag = Fragment.createFragment(gjf);

        Catalyst cat = new Catalyst(frag);
        System.out.println("Catalyst:\n" + cat);

        GJFfragment gjf2 = new GJFfragment("test2.gjf");
        Fragment frag2 = Fragment.createFragment(gjf2);

        cat = cat.addLeft(frag2);
        System.out.println("\n\nCatalyst:\n" + cat);

        GJFfragment gjf3 = new GJFfragment("test2.gjf");
        Fragment frag3 = Fragment.createFragment(gjf3);
        
        cat = cat.addRight(frag3);
        System.out.println("\n\nCatalyst:\n" + cat);
       
        MOL2InputFile file = new MOL2InputFile(cat);
        file.write("test_join.mol2");
    }
}
