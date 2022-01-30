import java.io.*;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

/**
 * @author hadziahmetovic on 23.01.22
 */
public class Reader {
    /*
     add GO reader
     add ensembl/gaf reader -> autoconvert/scrape info to generate ensembl verison
     diff-output reader -> flexible with col selector -> make gene, fc, label format??
     */



    private void readObo(){

        String line, id = null, name, namespace = null;
        Set<String> set = null;
        BufferedReader br = null;
        boolean good = true;

        try {br = new BufferedReader(new FileReader(obo));} catch (FileNotFoundException e) {}
        try {
            while ((line = br.readLine()) != null) {
                if (line.length() > 1) {
                    if (line.charAt(0) == '[' && line.charAt(2) == 'y') break;

                    else if (line.charAt(0) == 'i' && line.charAt(1) == 'd') {
                        if (!good) { }  //TODO nice if man
                        else if (set != null) {
                            goParent.put(id, set);
                            allGO.add(id);
                            allGO.addAll(set);
                        }
                        id = line.substring(4);
                        set = new HashSet<>();
                        good = true;
                    } else if (line.charAt(0) == 'n' && line.charAt(4) == ':') {
                        name = line.substring(6);
                        if (name.contains("obsolete")) good = false;
                        else goName.put(id, name);
                    } else if (line.charAt(0) == 'n' && line.charAt(4) == 's') {
                        namespace = line.substring(11);
                        if (!namespace.equals(root)) {
                            good = false;
                        }
                    } else if (line.charAt(1) == 's' && line.charAt(2) == '_' && line.charAt(3) == 'a') {
                        if (namespace.equals(root)) {
                            set.add(line.substring(6, line.indexOf("!", 7) - 1));
                        }
                    }
                }
            }
            br.close();

            if (set != null) {
                goParent.put(id, set);          //TODO this might be a problem because no check if good
            }
        } catch (IOException e) {
            System.out.println("IOexception at obo file");
        }
    }

    private void readMappringGAF(){

        int c = 0;
        String line;
        BufferedReader br = null;
        int one, two, three, four, five;
        String gene, go;

        try {br = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(mapping))));} catch (Exception e) {}
        try {
            while ((line = br.readLine()) != null){
                if (line.charAt(0) == '!') continue;
                one = line.indexOf("\t");
                two = line.indexOf("\t", one+1);
                three = line.indexOf("\t", two+1);
                four = line.indexOf("\t", three+1);
                five = line.indexOf("\t", four+1);
                //System.out.println(line);


                if (four-three > 1) {
                    //System.out.println(line);
                    continue;
                }



                gene = line.substring(two+1, three);

                if (geneSignif.get(gene) == null) {
                    continue;
                }

                go = line.substring(four+1, five);

                if (!allGO.contains(go)) {
                    continue;
                }

//                System.out.println("Gene\t"+gene);
//                System.out.println("Go\t"+go);
                //System.out.println(++c);
//                if (geneGO.get(line.substring(two+1, three)) == null) {
//                    Set<String> set = new HashSet<>();
//                    set.add(line.substring(four+1, five));
//                    geneGO.put(line.substring(two+1, three), set);
//                } else {
//                    geneGO.get(line.substring(two+1, three)).add(line.substring(four+1, five));
//                }

                if (!geneSignif.get(gene)) { //TODO my bad, should be true but its just for consistency now
                    if (goPositive.get(go) != null) {
                        goPositive.get(go).add(gene);
                    } else {
                        Set<String> set = new HashSet<>();
                        set.add(gene);
                        goPositive.put(go, set);
                    }
                } else {
                    if (goNegative.get(go) != null) {
                        goNegative.get(go).add(gene);
                    } else {
                        Set<String> set = new HashSet<>();
                        set.add(gene);
                        goNegative.put(go, set);
                    }
                }

            }
            br.close();
        }
        catch (IOException e) {
            System.out.println("IOexception gaf reader");
        }
    }

    private void readMappringEnsebl(){

        String line, name;
        BufferedReader br = null;
        int one, two;
        boolean b;

        try {br = new BufferedReader(new FileReader(mapping));} catch (FileNotFoundException e) {}
        try {
            line = br.readLine();
            while ((line = br.readLine()) != null) {
                one = line.indexOf("\t");
                two = line.indexOf("\t", one + 1);
                if (two - one < 2) continue;
                name = line.substring(one + 1, two);
                if (!geneSignif.keySet().contains(name)) {
                    continue;
                }
                //currentSet = new HashSet<>();
                currentGene = name;

                b = geneSignif.get(name);
                getGOs(line.substring(two + 1), b);
                //geneGO.put(name, currentSet);
            }
            br.close();
        }
        catch (IOException e) {
            System.out.println("IOexception ensembl reader");
        }
    }

    private void getGOs (String line, boolean b) {
        try {
            while (true) {
                if (!allGO.contains(line.substring(0, 10))) {
                    line = line.substring(11);
                    continue;
                }
                //currentSet.add(line.substring(0, 10));
                if (!b) {   //TODO dont ask why, it works
                    if (goPositive.get(line.substring(0, 10)) == null) {
                        Set<String> set123 = new HashSet<>();
                        set123.add(currentGene);
                        goPositive.put(line.substring(0, 10), set123);
                    } else {
                        goPositive.get(line.substring(0, 10)).add(currentGene);
                    }
                } else {
                    if (goNegative.get(line.substring(0, 10)) == null) {
                        Set<String> set123 = new HashSet<>();
                        set123.add(currentGene);
                        goNegative.put(line.substring(0, 10), set123);
                    } else {
                        goNegative.get(line.substring(0, 10)).add(currentGene);
                    }
                }
                line = line.substring(11);
            }
        } catch (Exception e) {}
    }

    private void readEnrich(){

        String line, name;
        BufferedReader br = null;
        int one, two;

        try {br = new BufferedReader(new FileReader(enrich));} catch (FileNotFoundException e) {}
        try {
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) == '#') isTrue.add(line.substring(1));
                else break;
            }

            while ((line = br.readLine()) != null) {
                one = line.indexOf("\t");
                two = line.indexOf("\t", one + 1);
                name = line.substring(0, one);
                if (geneFC.keySet().contains(name)) {
                    System.out.println(name);
                }
                geneFC.put(name, Double.parseDouble(line.substring(one + 1, two)));
                if (line.charAt(two + 1) == 'f') {
                    geneSignif.put(name, false);
                } else {
                    geneSignif.put(name, true);
                }
            }
            br.close();
        }
        catch (IOException e) {
            System.out.println("IOexception in enrich");
        }
    }

}
