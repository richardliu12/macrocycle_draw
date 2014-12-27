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

    /** 
     * Should be used to list items in the DATABASE.  Using DATABASE.toString
     * generates an excessively verbose output.
     */
     public static String getDatabase()
     {
         FragmentType[] types = {FragmentType.UREA, FragmentType.THIOUREA,
                                FragmentType.LINKER_1, FragmentType.LINKER_2,
                                FragmentType.LINKER_3, FragmentType.LINKER_4};
         String returnString = "";
         for ( FragmentType t : types )
            {
                returnString += "\n\n" + t + ": \n";
                for ( Fragment f : DATABASE.get(t))
                {
                    returnString += f.name + "\n";
                }
            }

         return returnString;
     }

    /**
    * Returns a List of cyclized Catalysts from the DATABASE. 
    * This requires a template of how to construct the catalyst,
    * given as a List of FragmentTypes.  These will have been
    * cyclized, but not yet minimized.  This method constructs
    * from LEFT TO RIGHT.
    * @param template
    * @return the list of Catalysts
    */
    public static List<Catalyst> createCatalysts(List<FragmentType> template)
    {
        // This will work by going through the template file,
        // and adding the appropriate fragments at each stage.
        // For each iteration, we populate the List nextCatalysts
        // based on currentCatalysts, and then move the newly
        // generated catalysts to currentCatalysts. 
        List<Catalyst> currentCatalysts = new ArrayList<>();
        List<Catalyst> nextCatalysts = new ArrayList<>();
          
        // Add the first fragment separately
        if ( template.size() == 0 )
            throw new IllegalArgumentException("Empty template!");
        else
            for ( Fragment f : DATABASE.get(template.get(0)) )
                currentCatalysts.add(new Catalyst(f));

        // Check that the catalyst list is not empty!
        if ( currentCatalysts.size() == 0 )
            throw new IllegalArgumentException("DATABASE contains no " + template.get(0));

        for ( int i = 1 ; i < template.size(); i++ )
            {
                for ( Catalyst c : currentCatalysts )
                    for ( Fragment f : DATABASE.get(template.get(i)) )
                        nextCatalysts.add(c.addRight(f));
                currentCatalysts = nextCatalysts;
                nextCatalysts = new ArrayList<>();
            }      

        // Cyclize the catalysts!
        for ( Catalyst c : currentCatalysts )
            nextCatalysts.add(c.cyclize());

        return nextCatalysts;
    }
}


