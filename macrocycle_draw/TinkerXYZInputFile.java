import java.io.*;
import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/** Class that creates xyz files from molecules. Has static method to convert Molecule class to xyz file **/
public class TinkerXYZInputFile extends InputFileFormat {
    
    /**
     * Constructor that writes a Molecule to a XYZ File.
     * @param molecule the input Molecule that will be converted to an xyz file
     */
    public TinkerXYZInputFile(Molecule molecule)
    {
	    super(molecule.toXYZString());
    }

}
