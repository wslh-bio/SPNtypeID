#!/usr/bin/env groovy

import java.util.zip.GZIPInputStream

// Function that actually looks at the file to see if it is empty or not
static def checkFastqFile(filePath) {
    File inputfile = new File(filePath)

    try {
        GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(inputfile))
        BufferedReader reader = new BufferedReader(new InputStreamReader(gzipInputStream))

        int lineCount = 0

        while ((reader.readLine()) != null) {
            lineCount++
            if (lineCount >= 4) {
                reader.close()
                return true
            }
        }
        reader.close()
    } catch (IOException e) {
        println "Error reading file: ${e.message}"
        return false
    }

    return false
}

// Starting function to process the samplesheet filepaths
static def processMeta( filePaths ) {

    if (!(filePaths instanceof List)) {
        filePaths = filePaths.toList()
    }

    // Assigns true or false based on the checkFastqFile function
    boolean pass = filePaths.every { filePath ->
        checkFastqFile(filePath as String) // Check each file and ensure all pass
    }

    return [pass]
}

// Function that writes the failed files to their own output file
static def failedFile(meta, outdir) {

    def fileName = "$outdir/Empty_samples.txt"
    def newFile = new File(fileName)

    // Extracts everything between the : and , to get the sample name
    def sample_name = (meta =~ /:(.+?),/).find() ? (meta =~ /:(.+?),/)[0][1] : null

    if (!newFile.exists()) {
        newFile.createNewFile()
    }

    newFile.withWriterAppend { writer ->
        writer.write("$sample_name\n")
    }
}
