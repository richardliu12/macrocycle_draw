import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
* Holds and organizes all fragments.  All the pieces in a Catalyst are
* called Fragments, and this Singleton reads all of those from the
* input directory.  The entries are stored in a map called DATABASE,
* and are sorted by fragment type.  Explanation courtesy of Prof. Kwan.
*
*/
public final class FragmentLibrary implements Singleton
{
    /** Do not try to instantiate this! */
    private FragmentLibrary() {throw new IllegalArgumentException("Not instantiable!");}
   
    /** Fragment collection, sorted by type. */
    public static final ImmutableMap<FragmentType,List<Fragment>> DATABASE;
    
    static{
        // collections of fragments by type
        List<Fragment> ureaList = new ArrayList<>();
        List<Fragment> thioureaList = new ArrayList<>();
        List<Fragment> linkerList_1 = new ArrayList<>();
        List<Fragment> linkerList_2 = new ArrayList<>();
        List<Fragment> linkerList_3 = new ArrayList<>();
        List<Fragment> linkerList_4 = new ArrayList<>();

        File fragmentDirectory = new File(Settings.INPUT_DIRECTORY);
        for ( File f:fragmentDirectory.listFiles() )
        {
            String filename = f.getName();
            if ( filename.endsWith(".gjf") )
            {
                GJFfragment gjf = new GJFfragment(Settings.INPUT_DIRECTORY+filename);
                Fragment fragment = Fragment.createFragment(gjf);
                switch(fragment.fragmentType)
                {
                    case UREA: ureaList.add(fragment);
                        break;
                    case THIOUREA: thioureaList.add(fragment);
                        break;
                    case LINKER_1: linkerList_1.add(fragment);
                        break;
                    case LINKER_2: linkerList_2.add(fragment);
                        break;
                    case LINKER_3: linkerList_3.add(fragment);
                        break;
                    case LINKER_4: linkerList_4.add(fragment);
                        break;
                    default:
                        throw new IllegalArgumentException("Fragment type error in reading " + filename);
                }
            }
        }
        
        Map<FragmentType,List<Fragment>> tempMap = new HashMap<>();
        tempMap.put(FragmentType.UREA, ureaList);
        tempMap.put(FragmentType.THIOUREA, thioureaList);
        tempMap.put(FragmentType.LINKER_1, linkerList_1);
        tempMap.put(FragmentType.LINKER_2, linkerList_2);
        tempMap.put(FragmentType.LINKER_3, linkerList_3);
        tempMap.put(FragmentType.LINKER_4, linkerList_4);

        DATABASE = ImmutableMap.copyOf(tempMap);
    }
}


