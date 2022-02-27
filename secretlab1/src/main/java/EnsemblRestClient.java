import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.net.HttpURLConnection;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.Reader;
import net.minidev.json.*;
import net.minidev.json.parser.*;

public class EnsemblRestClient {

    public static final String SERVER = "http://rest.ensembl.org";
    public static final JSONParser PARSER = new JSONParser(JSONParser.MODE_JSON_SIMPLE);

    public static int requestCount = 0;
    public static long lastRequestTime = System.currentTimeMillis();


    public static void main(String[] args) throws Exception {

        System.out.println(">ENSG00000105939");
        System.out.println(getSequence("ENSG00000105939"));

        System.out.println(getGeneID("homo_sapiens", "ZC3HAV1"));

        getRegions("ENSP00000242351", 1, 902);
    }

    public static void testVariants() throws ParseException, InterruptedException, IOException {
        String species, symbol;
        species = "human";
        symbol = "BRAF";

        JSONArray variants = getVariants(species, symbol);
        for(Object variantObject: variants) {
            JSONObject variant = (JSONObject)variantObject;
            String srName = (String)variant.get("seq_region_name");
            Number start = (Number)variant.get("start");
            Number end = (Number)variant.get("end");
            Number strand = (Number)variant.get("strand");
            String id = (String)variant.get("id");
            String consequence = (String)variant.get("id");
            String output = String.format("%s:%d-%d:%d ==> %s (%s)", srName, start, end, strand, id, consequence);
            System.out.println(output);
        }
    }

    public static void getRegions(String protein_id, int start, int end) throws IOException {
        String server = "http://rest.ensembl.org";
        String ext = "/map/translation/"+protein_id+"/"+start+".."+end+"?";
        URL url = new URL(server + ext);

        URLConnection connection = url.openConnection();
        HttpURLConnection httpConnection = (HttpURLConnection) connection;

        httpConnection.setRequestProperty("Content-Type", "application/json");


        InputStream response = connection.getInputStream();
        int responseCode = httpConnection.getResponseCode();

        if (responseCode != 200) {
            throw new RuntimeException("Response code was not 200. Detected response was " + responseCode);
        }

        String output;
        Reader reader = null;
        try {
            reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
            StringBuilder builder = new StringBuilder();
            char[] buffer = new char[8192];
            int read;
            while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
                builder.append(buffer, 0, read);
            }
            output = builder.toString();
        } finally {
            if (reader != null) try {
                reader.close();
            } catch (IOException logOrIgnore) {
                logOrIgnore.printStackTrace();
            }
        }

        System.out.println(output);
    }


    /**
     * http://rest.ensembl.org/documentation/info/xref_external
     * @param species Species name/alias. Examples: homo_sapiens, human
     * @param symbol Symbol or display name of a gene. Example: BRCA2
     * @return gene_id matching the input symbol
     */
    public static String getGeneID(String species, String symbol) throws ParseException, IOException, InterruptedException {
        String endpoint = "/xrefs/symbol/"+species+"/"+symbol+"?object_type=gene";
        JSONArray genes = (JSONArray) getJSON(endpoint);
        if(genes.isEmpty()) {
            throw new RuntimeException("Got nothing for endpoint "+endpoint);
        }
        JSONObject gene = (JSONObject)genes.get(0);
        return (String)gene.get("id");
    }

    /**
     * Converts the content of a request into JSON
     * @param endpoint GET request. Example: /xrefs/symbol/homo_sapiens/BRCA2?object_type=gene
     * @return Object which can be cast to JSONObject in order to access content.
     */
    public static Object getJSON(String endpoint) throws ParseException, IOException, InterruptedException {
        String jsonString = getContent(endpoint);
        return PARSER.parse(jsonString);
    }

    public static String getSequence(String id) throws ParseException, InterruptedException, IOException {
        String endpoint = "/sequence/id/" +id + "?";
        JSONObject entry = (JSONObject) getJSON(endpoint);
        return (String) entry.get("seq");
    }

    public static JSONArray getVariants(String species, String symbol) throws ParseException, IOException, InterruptedException {
        String id = getGeneID(species, symbol);
        return (JSONArray) getJSON("/overlap/id/"+id+"?feature=variation");
    }

    public static String getContent(String endpoint) throws IOException, InterruptedException {

        if(requestCount == 15) { // check every 15
            long currentTime = System.currentTimeMillis();
            long diff = currentTime - lastRequestTime;
            //if less than a second then sleep for the remainder of the second
            if(diff < 1000) {
                Thread.sleep(1000 - diff);
            }
            //reset
            lastRequestTime = System.currentTimeMillis();
            requestCount = 0;
        }

        URL url = new URL(SERVER+endpoint);
        URLConnection connection = url.openConnection();
        HttpURLConnection httpConnection = (HttpURLConnection)connection;
        httpConnection.setRequestProperty("Content-Type", "application/json");
        requestCount++;

        InputStream response = httpConnection.getInputStream();
        int responseCode = httpConnection.getResponseCode();

        if(responseCode != 200) {
            if(responseCode == 429 && httpConnection.getHeaderField("Retry-After") != null) {
                double sleepFloatingPoint = Double.valueOf(httpConnection.getHeaderField("Retry-After"));
                double sleepMillis = 1000 * sleepFloatingPoint;
                Thread.sleep((long)sleepMillis);
                return getContent(endpoint);
            }
            throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
        }

        String output;
        Reader reader = null;
        try {
            reader = new BufferedReader(new InputStreamReader(response, "UTF-8"));
            StringBuilder builder = new StringBuilder();
            char[] buffer = new char[8192];
            int read;
            while ((read = reader.read(buffer, 0, buffer.length)) > 0) {
                builder.append(buffer, 0, read);
            }
            output = builder.toString();
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException logOrIgnore) {
                    logOrIgnore.printStackTrace();
                }
            }
        }

        return output;
    }
}
