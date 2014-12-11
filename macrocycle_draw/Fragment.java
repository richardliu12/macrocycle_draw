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
}
