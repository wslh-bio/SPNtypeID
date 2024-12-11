import java.util.zip.GZIPInputStream

class CheckGZIPLines {

    // Create a true or false scenario where true is data and false is no data
    // Takes in the file paths as a string
    static boolean DataInTheFirstFourLines(String filePath) {
        // First try this chunk
        try {

            // Sets a file variable to the file path string that was entered into the class 
            File file = new File(filePath)

            // If the file does not exist OR the file's length = 0
            if (!file.exists() || file.length() == 0) {

                // returns a false 
                return false // File is missing or empty
            }

            // new FileInputStream creates a stream of RAW bytes from the specific file, these bytes need to be decompressed
            // new GZIPInputStream(new FileInputStream(file)) will actually decompress the file as you read it in, decompresses the raw bytes into bytes
            // Stores gzipInputStream as a variable with type GZIPInputStream
            GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(file))

            // new InputStreamReader(gzipInputStream) reads the bytes and then makes them into text
            // new BufferedReader(new InputStreamReader(gzipInputStream)) goes line by line 
            // Stores reader as a variable with type BufferedReader
            BufferedReader reader = new BufferedReader(new InputStreamReader(gzipInputStream))

            // Sets up an interger variable called lineCount to 0
            int lineCount = 0

            // Creates a new string variable which is line
            String line

            // While line(which is equal to reading the line in the file) does not = nothing and is less than 4
            while ((line = reader.readLine()) != null && lineCount < 4) {
                // If not the trimmed line is empty 
                if (!line.trim().isEmpty()) {
                    // stop reading the file
                    reader.close()
                    // Found a non-empty line within the first 4 lines
                    return true 
                }
                // Adds 1 to line count variable so you can iterate through the while loop
                lineCount++
            }

            // Closes the reader and returns false if the first 4 lines are 
            reader.close()
            return false // All first 4 lines are empty

        // Allows to catch this exception
        } catch (Exception e) {
            println("Error reading file ${filePath}: ${e.message}")
            return false
        }
    }
}