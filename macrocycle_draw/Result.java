import java.io.*;

public interface Result extends Serializable
{
    public static final long serialVersionUID = 1L;

    public static final Result JOB_FAILED = new Result() {
        public static final long serialVersionUID = 1L;
        public String toString() { return "job failed"; } };

    public static final Result JOB_INTERRUPTED = new Result() {
        public static final long serialVersionUID = 1L;
        public String toString() { return "job was interrupted"; } };

    public static final Result JOB_UNAVAILABLE = new Result() {
        public static final long serialVersionUID = 1L;
        public String toString() { return "job result not available"; } };

    public static final Result JOB_COMPLETE = new Result() {
        public static final long serialVersionUID = 1L;
        public String toString() { return "job complete"; } };

    public String toString();
}
