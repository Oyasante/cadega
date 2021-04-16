/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package padll;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import static padll.Main.DRUGNAME;
import static padll.Main.DRUGS;

/**
 *
 * @author Yaw's PC
 */
public class Integrator {
    public static final String PATH = "";
    public static final String PATH3 = "..\\predata\\";
    public static final String PATH5 = "..\\drugs\\";
    public static final String PATH6 = "";
    
    
    public static void main (String[] args) throws IOException {
        scanDrugSet();
        //stripNetwork();
        
        String[] GENES_NP1 = {"TP53_mut","PPM1D_amp","CSNK2A1_del", "SLC25A11_mut"};
        String[] GENES_CD1 = {"SLC25A11_mut"};
   
        String[] GENES_NP2 = {"CCND1_mut","PSMD3_amp","CDKN2A_del","HCCS_mut"};
        String[] GENES_CD2 = {"HCCS_mut"};
        
        String[] GENES_NP3 = {"TNS3_mut","ARHGEF4_mut","DLC1_del","NPPA_mut"};
        String[] GENES_CD3 = {"NPPA_mut"};
    
        String[] GENES_NP4 = {"CTNNB1_mut","SOX9_mut","AR_del","YIF1B_amp"};
        String[] GENES_CD4 = {"YIF1B_amp"};
    
        String[] GENES_NP5 = {"DNAH7_mut","TPR_mut","UPF3B_del","ZBTB8OS_mut"};
        String[] GENES_CD5 = {"ZBTB8OS_mut"};
    
        String[] GENES_NP6 = {"KRAS_mut","CCND1_amp","RB1_del","ZNF385A_del"};
        String[] GENES_CD6 = {"ZNF385A_del"};
        
        selectNetwork("NPDrug1", GENES_NP1);
        //selectNetwork("CDDrug1", GENES_CD1);
        
        selectNetwork("NPDrug2", GENES_NP2);
        //selectNetwork("CDDrug1", GENES_CD2);
        
        selectNetwork("NPDrug3", GENES_NP3);
        //selectNetwork("CDDrug1", GENES_CD3);
        
        selectNetwork("NPDrug4", GENES_NP4);
        //selectNetwork("CDDrug1", GENES_CD4);
        
        selectNetwork("NPDrug5", GENES_NP5);
        //selectNetwork("CDDrug1", GENES_CD5);
        
        selectNetwork("NPDrug6", GENES_NP6);
        //selectNetwork("CDDrug1", GENES_CD6);
    }
    public static void scanDrugSet() throws FileNotFoundException, IOException {
        ArrayList<ArrayList<TissueLink>> allVals = new ArrayList<>();
        ArrayList<String> maxTis = new ArrayList<>();
        
        for (int i = 0; i < DRUGS.length; i++) {
            File f = new File(PATH3+ DRUGS[i] + ".csv");
            Scanner sc = new Scanner(f);
            sc.useDelimiter("\n");
            int outerLoopCount = 0;
            ArrayList<TissueLink> tis = new ArrayList<>();
            //ArrayList<String> curTissues = new ArrayList<>();

            while (sc.hasNext()) {
                Scanner sc2 = new Scanner(sc.next());
                sc2.useDelimiter(",");
                int innerLoopCount = 0;
                String cellline = "";
                TissueLink cur = new TissueLink();
                

                while (sc2.hasNext()) {
                    String current = sc2.next();


                    if (innerLoopCount == 0 && outerLoopCount == 0) {
                        innerLoopCount++;
                        continue;
                    } else if (innerLoopCount == 0 && outerLoopCount > 0){
                        cellline = current;
                    } else if (innerLoopCount > 0 && outerLoopCount == 0) {
                        innerLoopCount++;
                        continue;
                    } else if (innerLoopCount > 0 && outerLoopCount > 0) {
                        
                        if (innerLoopCount == 3) {
                            cellline = processCellLine(cellline);
                            String totalTissue = cellline + "_" + current;
                            cur.setTissue(totalTissue);
                            cur.setDrug(DRUGNAME[i]);
                            //curTissues.add(totalTissue);
                            if (!(maxTis.contains(totalTissue)))
                                maxTis.add(totalTissue);

                        }
                        if (innerLoopCount == 5) {
                            double value = Double.parseDouble(current);
                            cur.setValue(value);
                        }

                    }
                    innerLoopCount++;
                }
                tis.add(cur);
                
                outerLoopCount++;
            }
            System.out.println("Done with line " + i + "!");
            allVals.add(tis);
        }
        
        StringBuilder header = new StringBuilder("Description\t");
        boolean first = true;
        
        for (int j = 0; j < DRUGS.length; j++) {
            ArrayList<TissueLink> active = allVals.get(j);
            String drug = active.get(0).getDrug();
            
            StringBuilder row = new StringBuilder("");
            
            // list of values to update
            ArrayList<Double> toBeNormed = new ArrayList<>();
            double mean = 0;
            int count = 0;
            
            for (String cur : maxTis) {
                if (first) header.append(cur.toUpperCase()).append("\t");
                
                int index = findTis(active, cur);
                
                if (index > 0) {
                    double val = active.get(index).getSensitivity();
                    toBeNormed.add(val);
                    if (!Double.isNaN(val))
                        mean += val;
                    row.append("" + val).append("\t");
                } else {
                    toBeNormed.add(Double.NaN);
                    row.append("NaN").append("\t");
                }
                count++;
            }
            if (first)
                first = false;
            
            File f1 = new File(PATH5 + "NPTarget" + j + ".txt");
            FileWriter fw = new FileWriter(f1);
            fw.write(header.toString());
            fw.write("\n");
            fw.write(DRUGNAME[j] + "\t");
            fw.write(row.toString());
            fw.flush();
            fw.close();
            
            // run Z norm
            mean = mean/count;
            ArrayList<Double> zScores = zNorm(toBeNormed, mean);
            
            File f2 = new File(PATH5 + "CDTarget" + j + ".txt");
            FileWriter fw2 = new FileWriter(f2);
            fw2.write(header.toString());
            fw2.write("\n");
            fw2.write(DRUGNAME[j] + "\t");
            fw2.write(getRow(zScores));
            fw2.flush();
            fw2.close();
            
        }
        
    }
    
    public static int findTis(ArrayList<TissueLink> a, String c) {
        int ret = -1;
        int count = 0;
        for (TissueLink t: a) {
            if (c.equals(t.getTissue())){
                ret = count;
                return ret;
            }
            count++;
        }
        return ret;
    }

    private static String processCellLine(String cellline) {
        String ret = "";
        for (int i = 0; i < cellline.length();i++) {
            char t = cellline.charAt(i);
            if( t != '-') {
                ret += t;
            }
        }
        return ret;
    }

    // one-time operation on STRING Database, do not use without changing PATH
    private static void stripNetwork() throws FileNotFoundException, IOException {
        File f = new File(PATH + "HumanStringNet.txt");
            Scanner sc = new Scanner(f);
            sc.useDelimiter("\n");
            int outerLoopCount = 0;
            StringBuilder network = new StringBuilder();
            //ArrayList<String> curTissues = new ArrayList<>();

            while (sc.hasNext()) {
                Scanner sc2 = new Scanner(sc.next());
                sc2.useDelimiter("\t");
                int innerLoopCount = 0;
                
                
                while (sc2.hasNext()) {
                    if (innerLoopCount >= 2) break;
                    network.append(sc2.next());
                    if (innerLoopCount < 1)
                        network.append("\t");
                    innerLoopCount++;
                }
                network.append("\n");
            }
            File f2 = new File(PATH6 + "network.tsv");
            FileWriter fw = new FileWriter(f2);
            fw.write(network.toString());
            fw.flush();
            fw.close();
    }

    private static ArrayList<Double> zNorm(ArrayList<Double> toBeNormed, double mean) {
        ArrayList<Double> ret = new ArrayList<>();
        int len = toBeNormed.size();
        
        double sum = 0;
        for (int j = 0; j < len; j++) {
            if (!Double.isNaN(toBeNormed.get(j)))
                sum += (toBeNormed.get(j) - mean)*(toBeNormed.get(j) - mean);
        }
        double std = Math.sqrt(sum / (len - 1));
        
        for (int j = 0; j < len; j++) {
            ret.add((toBeNormed.get(j) - mean)/std);
        }
        
        return ret;
    }

    private static String getRow(ArrayList<Double> zScores) {
        StringBuilder ret = new StringBuilder();
        for (Double d: zScores)
               ret.append("" + d).append("\t");
        return ret.toString();
    }

    private static void selectNetwork(String s, String[] check) throws FileNotFoundException, IOException {
         File f = new File(PATH + "HumanStringNet.txt");
            Scanner sc = new Scanner(f);
            sc.useDelimiter("\n");
            int outerLoopCount = 0;
            StringBuilder network = new StringBuilder();
            //ArrayList<String> curTissues = new ArrayList<>();

            while (sc.hasNext()) {
                Scanner sc2 = new Scanner(sc.next());
                sc2.useDelimiter("\t");
                int innerLoopCount = 0;
                boolean adder = false;
                String temp = "";
                
                while (sc2.hasNext()) {
                    if (innerLoopCount >= 2) break;
                    String val = (sc2.next());
                    
                    for (int j = 0; j < check.length; j++) {
                        String b = check[j].substring(0, check[j].indexOf("_"));
                        if (b.equals(val)) {
                            adder = true;
                            break;
                        }
                            
                    }
                    temp = temp + val;
                    if (innerLoopCount < 1)
                        temp = temp + "\t";
                    innerLoopCount++;
                }
                if (adder)
                    network.append(temp + "\n");
            }
            File f2 = new File(PATH6 + "selectNet" + s + ".tsv");
            FileWriter fw = new FileWriter(f2);
            fw.write(network.toString());
            fw.flush();
            fw.close();
    }
}
