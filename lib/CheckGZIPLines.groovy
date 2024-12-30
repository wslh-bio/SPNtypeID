#!/usr/bin/env groovy
import java.util.zip.GZIPInputStream
import nextflow.Channel
import nextflow.Nextflow

class CheckGZIPLines {

    public static boolean DataInTheFirstFourLines(String filePath) {
        try {
            File file = new File(filePath)
            if (!file.exists() || file.length() == 0) {
                return false // File is missing or empty
            }
            GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(file))
            BufferedReader reader = new BufferedReader(new InputStreamReader(gzipInputStream))

            int lineCount = 0
            String line
            while ((line = reader.readLine()) != null && lineCount < 4) {
                if (!line.trim().isEmpty()) {
                    reader.close()
                    return true // Found a non-empty line
                }
                lineCount++
            }

            reader.close()
            return false // All first 4 lines are empty
        } catch (Exception e) {
            println("Error reading file ${filePath}: ${e.message}")
            return false
        }
    }
}