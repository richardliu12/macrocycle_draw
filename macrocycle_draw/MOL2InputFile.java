/** Represents a MOL2 file for debugging purposes. Designed to be PYMOL-readable.*/
public class MOL2InputFile extends InputFileFormat
{
    public MOL2InputFile(Molecule molecule)
    {
        super(molecule.toMOL2());
    }
}
