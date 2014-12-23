/**
* Describes the type of molecular fragment.  Urea and Thiourea 
* are self-explanatory, but there are several types of linkers.
* This is designed to allow for different types linkers, so that
* fragments in different linker groups will not mix.
*/

public enum FragmentType{
    UREA, THIOUREA, LINKER_1, LINKER_2, LINKER_3, LINKER_4
}
