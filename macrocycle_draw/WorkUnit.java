import java.io.*;
import java.util.concurrent.*;

public interface WorkUnit extends Callable<Result>
{
    public Result call();
    public String toString();
}
