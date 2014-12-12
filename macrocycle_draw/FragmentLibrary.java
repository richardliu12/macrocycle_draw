import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
* Holds and organizes all fragments.  Singleton class.
*
*/
public class FragmentLibrary implements Singleton
{
    private static FragmentLibrary INSTANCE = null;
    protected FragmentLibrary {};

    // collections of fragments by type
    public ImmutableList<Fragment> ureaList;
    public ImmutableList<Fragment> thioureaList;
    public ImmutableList<Fragment> linkerList_1;
    public ImmutableList<Fragment> linkerList_2;
    public ImmutableList<Fragment> linkerList_3;
    public ImmutableList<Fragment> linkerList_4;

    /** Creates lone instance of FragmentLibrary, may only be called once.
    */
    public static FragmentLibrary getInstance()
    {
        if ( INSTANCE == null )
            INSTANCE = new FragmentLibrary;
        else
            throw new IllegalArgumentException("FragmentLibrary is a singleton!");
        return INSTANCE;
    }

    public void populate(String[] filenames)
    {
    }
}


