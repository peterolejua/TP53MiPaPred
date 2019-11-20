
/**
 *
 * Copyright (c) @author Xiaoming Liu, Ph.D.
 * Associate Professor,
 * USF Genomics,
 * College of Public Health,
 * University of South Florida at Tampa
 * 
 * This source code is distributed under the RECEX SHARED SOURCE LICENSE
 * 
 * You are free to download, copy, compile, study, and refer to the source code for any personal use of yours.
 * You are free to make any modifications to the source covered by this license.
 * You may NOT under any circumstance copy, redistribute and/or republish the source or a work based on it (which
 * includes binary or object code compiled from it) in part or whole.
 * If you intend to incorporate the source code, in part or whole, into any free or proprietary program, you need to explicitly
 * write to the original author(s) to ask for permission.
 * The source code licensed under this license is shared "as is".
 * 
 * You shall already get a copy of the license, if not, you can obtain a copy at 
 * https://raw.github.com/Recex/Licenses/master/SharedSourceLicense/LICENSE.txt
 */
import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
public class search_dbNSFP40b1a {
    public static void main(String[] args)throws Exception{
        String infile="";
        String outfile="";
        String outfile2="";
        String classdir=(new File(search_dbNSFP40b1a.class.getProtectionDomain().getCodeSource().getLocation().getPath())).getCanonicalPath()+File.separator;
        boolean hg18=false;
        boolean hg19=false;
        boolean hg38=true;
        boolean allchr=true;
        boolean searchallchr=false;
        boolean hasgene=false;
        boolean hasens=false;
        boolean hasuniprot=false;
        boolean hasenspsnp=false;
        boolean hassnp=false;
        boolean hasentrez=false;
        boolean hasuniprotsnp=false;
        boolean hasrs=false;
        boolean allcol=true;
        boolean isvcf=false;
        int[] cols=new int[0];
        boolean preservevcf=false;
        boolean searchsrsnp=false;
        HashMap foundgene=new HashMap();
        HashMap unfoundgene=new HashMap();
        String[] hg18chr={"6","11","14","5","8"};
        String[] hg19chr={"Y","11","14","5","4"};
        String[] hg38chr={"6","17","22","1","4"};
        String[] chrs={"M","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};
        HashMap chrsx=new HashMap();
        if(args.length==0||args.length<4){
            System.out.println("Usage: search_dbNSFP40b1a [Options] -i input_file_name -o output_file_name");
            System.out.println("Options (not necessary):");
            System.out.println("-v human_genome_version");
            System.out.println("   human_genome_version=hg18, hg19 or hg38, default=hg38");
            System.out.println("   for example \"-v hg19\"");
            System.out.println("-c list_of_chromosomes_to_search");
            System.out.println("   chromosomes are represented by 1,2,...,22,X,Y, separated by commas");
            System.out.println("   default=all chromosomes");
            System.out.println("   for example \"-c 1,3,22,X\"");
            System.out.println("-w list_of_columns_to_write_to_output_file");
            System.out.println("   columns are represented by 1,2,..., separated by commas,");
            System.out.println("   continuous number block can be simplified by begin-end");
            System.out.println("   default=all columns");
            System.out.println("   for example \"-w 1-6,8\"");
            System.out.println("-p");
            System.out.println("   preserve all columns of the input vcf file in the output file");
            System.out.println("   only work if the input file is in vcf format");
            System.out.println("   may require large memory for a large vcf file");
            System.out.println("-s");
            System.out.println("   will search attached database(s) dbscSNV (and SPIDEX, if available)");
            System.out.println("   will be outputed to output_file_name.dbscSNV (and output_file_name.SPIDEX)");
            System.exit(0);
        }
        else{
            for(int i=0;i<args.length;i++){
                if(args[i].equalsIgnoreCase("-v")&&args[i+1].equalsIgnoreCase("hg18")) {
                    hg18=true;
                    hg38=false;
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-v")&&args[i+1].equalsIgnoreCase("hg19")) {
                    hg19=true;
                    hg38=false;
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-v")&&args[i+1].equalsIgnoreCase("hg38")) {
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-c")){
                    StringTokenizer t=new StringTokenizer(args[i+1],",");
                    int nt=t.countTokens();
                    allchr=false;
                    for(int j=0;j<nt;j++){
                        String chr=t.nextToken();
                        chrsx.put(chr,"");
                    }
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-i")) {
                    infile=args[i+1];
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-o")) {
                    outfile=args[i+1];
                    i++;
                    if(searchsrsnp)outfile2=outfile;
                }
                else if(args[i].equalsIgnoreCase("-w")){
                    StringTokenizer t=new StringTokenizer(args[i+1],",");
                    int nt=t.countTokens();
                    allcol=false;
                    int nc=0;
                    for(int j=0;j<nt;j++){
                        String tmp=t.nextToken();
                        StringTokenizer t2=new StringTokenizer(tmp,"-");
                        if(t2.countTokens()==1)nc++;
                        else{
                            int begin=Integer.parseInt(t2.nextToken());
                            int end=Integer.parseInt(t2.nextToken());
                            nc+=end-begin+1;
                        }
                    }
                    cols=new int[nc];
                    t=new StringTokenizer(args[i+1],",");
                    nt=t.countTokens();
                    nc=0;
                    for(int j=0;j<nt;j++){
                        String tmp=t.nextToken();
                        StringTokenizer t2=new StringTokenizer(tmp,"-");
                        if(t2.countTokens()==1){
                            cols[nc]=Integer.parseInt(t2.nextToken());
                            nc++;
                        }
                        else{
                            int begin=Integer.parseInt(t2.nextToken());
                            int end=Integer.parseInt(t2.nextToken());
                            for(int k=begin;k<=end;k++){
                                cols[nc]=k;
                                nc++;
                            }
                        }
                    }
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-s")) {
                    searchsrsnp=true;
                    outfile2=outfile;
                    
                }
                else if(args[i].equalsIgnoreCase("-p")) {
                    isvcf=true;
                    preservevcf=true;
                }
                else {
                    System.out.println("wrong options!");
                    System.exit(0);
                }
            }
        }
        PrintWriter out=new PrintWriter(new FileWriter(outfile));
        PrintWriter out2=new PrintWriter(new FileWriter(outfile+".err"));

        //read input file
        HashMap snps4=new HashMap();
        HashMap snps6=new HashMap();
        HashMap snpsp=new HashMap();
        HashMap snprs=new HashMap();
        HashMap genes=new HashMap();
        HashMap poss=new HashMap();
        HashMap genequery=new HashMap();
        HashMap uniprotsnp=new HashMap();
        HashMap enspsnp=new HashMap();
        BufferedReader in=new BufferedReader(new FileReader(infile));
        if(infile.substring(infile.length()-4).equalsIgnoreCase(".vcf"))isvcf=true;
        if(infile.toLowerCase().endsWith(".gz")){
            in=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile))));
            if(infile.substring(infile.length()-7).toLowerCase().equals(".vcf.gz"))isvcf=true;
        }
        String title="";
        while(in.ready()){
            String line=in.readLine();
            if(line.isEmpty()) continue;
            if(line.charAt(0)=='#'){
                title=line;
                continue;
            }
            StringTokenizer t=new StringTokenizer(line);
            int nt=t.countTokens();
            if(isvcf){
                t=new StringTokenizer(line,"\t");
                String chr=t.nextToken().toUpperCase();
                if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("CHR"))chr=chr.substring(3);
                if(chr.equals("MT"))chr="M";
                if(chr.equals("XY"))chr="X";
                String pos=t.nextToken();
                String id=t.nextToken();
                String ref=t.nextToken().toUpperCase();
                String alt=t.nextToken().toUpperCase();
                if(preservevcf)snps4.put(chr+"_"+pos+"_"+ref+"_"+alt, compress(" "+line.trim()));
                else snps4.put(chr+"_"+pos+"_"+ref+"_"+alt, " ");
                if(allchr)chrsx.put(chr,"");
                hassnp=true;
            }
            else if(nt==6){
                String chr=t.nextToken().toUpperCase();
                if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("CHR"))chr=chr.substring(3);
                if(chr.equals("MT"))chr="M";
                if(chr.equals("XY"))chr="X";
                String pos=t.nextToken();
                String ref=t.nextToken().toUpperCase();
                String alt=t.nextToken().toUpperCase();
                String aaref=t.nextToken().toUpperCase();
                String aaalt=t.nextToken().toUpperCase();
                snps6.put(chr+"_"+pos+"_"+ref+"_"+alt+"_"+aaref+"_"+aaalt, " ");
                if(searchsrsnp)snps4.put(chr+"_"+pos+"_"+ref+"_"+alt, " ");
                if(allchr)chrsx.put(chr,"");
                hassnp=true;
            }
            else if(nt == 4)
            {
                String chr=t.nextToken().toUpperCase();
                if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("CHR"))chr=chr.substring(3);
                if(chr.equals("MT"))chr="M";
                if(chr.equals("XY"))chr="X";
                String pos=t.nextToken();
                String ref=t.nextToken().toUpperCase();
                String alt=t.nextToken().toUpperCase();
                snps4.put(chr+"_"+pos+"_"+ref+"_"+alt, " ");
                if(allchr)chrsx.put(chr,"");
                hassnp=true;
            }
            else if(nt==1){
                String tmp=t.nextToken().toUpperCase();
                if(tmp.startsWith("RS")&&isInt(tmp.substring(2))){
                    snprs.put(tmp.toLowerCase()," ");
                    //hassnp=true;
                    hasrs=true;
                }
                else if(tmp.length()>7&&tmp.substring(0,7).equalsIgnoreCase("Uniprot")){
                    StringTokenizer t2=new StringTokenizer(tmp,":");
                    int nt2=t2.countTokens();
                    if(nt2==2){
                        genequery.put(tmp,tmp);
                        hasgene=true;
                    }
                    else if(nt2==3){
                        snpsp.put(tmp," ");
                        StringTokenizer t3=new StringTokenizer(tmp,":");
                        t3.nextToken();
                        uniprotsnp.put(t3.nextToken(),"");//put id
                        hassnp=true;
                        hasuniprotsnp=true;
                    }
                    hasuniprot=true;
                }
                else if(tmp.length()>7&&tmp.substring(0,7).equalsIgnoreCase("Ensembl")){
                    StringTokenizer t2=new StringTokenizer(tmp,":");
                    int nt2=t2.countTokens();
                    if(nt2==2){
                        genequery.put(tmp,tmp);
                        hasgene=true;
                        hasens=true;
                    }
                    else if(nt2==3){
                        snpsp.put(tmp," ");
                        StringTokenizer t3=new StringTokenizer(tmp,":");
                        t3.nextToken();
                        enspsnp.put(t3.nextToken(),"");//put id
                        hassnp=true;
                        hasenspsnp=true;
                    }
                    
                }
                else if(tmp.length()>6&&tmp.substring(0,6).equalsIgnoreCase("Entrez")){
                    hasentrez=true;
                    genequery.put(tmp, tmp);
                    hasgene=true;
                }
                else {
                    genequery.put(tmp, tmp);
                    hasgene=true;
                }
            }
            else if(nt==2){
                String chr=t.nextToken().toUpperCase();
                if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("CHR"))chr=chr.substring(3);
                if(chr.equals("MT"))chr="M";
                if(chr.equals("XY"))chr="X";
                String pos=t.nextToken();
                poss.put(chr+"_"+pos, "");
                if(allchr)chrsx.put(chr,"");
                hassnp=true;
            }
        }
        in.close();
        if(hg19){//some snp in hg19 are mapped to a different chr in hg38
            for(int i=0;i<hg19chr.length;i++){
                if(chrsx.containsKey(hg19chr[i]))chrsx.put(hg38chr[i],"");
            }
        }
        else if(hg18){//some snp in hg18 are mapped to a different chr in hg38
            for(int i=0;i<hg18chr.length;i++){
                if(chrsx.containsKey(hg18chr[i]))chrsx.put(hg38chr[i],"");
            }
        }
        //search dbNSFP_gene
        int ngenecol=0;
        HashMap dbNSFPgene=new HashMap(33000);
        String genetitle="";
        String infile2=classdir+"dbNSFP4.0b1_gene.gz";
        File f = new File(infile2);
        if (!f.isFile())System.err.println("error: cannot find dbNSFP4.0b1_gene.gz");
        else{
            in=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile2))));
            String line=in.readLine();
            int indx=line.indexOf("chr");
            genetitle=line.substring(indx+4);
            while(in.ready()){
                line=in.readLine();
                StringTokenizer t=new StringTokenizer(line,"\t");
                int nt=t.countTokens();
                ngenecol=nt;
                String[] work=new String[nt];
                for(int i=0;i<nt;i++)work[i]=t.nextToken();
                String genename=work[0].toUpperCase();
                String chr=work[2];
                if(hasgene){
                    if(genequery.containsKey(genename)){
                        if(allchr)chrsx.put(chr,"");
                        genes.put(genename,genename);
                        genequery.remove(genename);
                    }
                    
                }
                if(hasens){
                    t=new StringTokenizer(work[1],";");
                    nt=t.countTokens();
                    for(int i=0;i<nt;i++){
                        String tmp=t.nextToken().toUpperCase();
                        if(genequery.containsKey("ENSEMBL:"+tmp)){
                            if(allchr)chrsx.put(chr,"");
                            genes.put("ENSEMBL:"+tmp,"ENSEMBL:"+tmp);
                            genes.put(genename,"ENSEMBL:"+tmp);
                            genequery.remove("ENSEMBL:"+tmp);
                        }
                    }
                }
                if(hasuniprot){
                    t=new StringTokenizer(work[5],";");
                    nt=t.countTokens();
                    for(int i=0;i<nt;i++){
                        String tmp=t.nextToken().toUpperCase();
                        if(genequery.containsKey("UNIPROT:"+tmp)){
                            if(allchr)chrsx.put(chr,"");
                            genes.put("UNIPROT:"+tmp,"UNIPROT:"+tmp);
                            genes.put(genename,"UNIPROT:"+tmp);
                            genequery.remove("UNIPROT:"+tmp);
                        }
                        if(hasuniprotsnp){
                            if(uniprotsnp.containsKey(tmp)){
                                if(allchr)chrsx.put(chr,"");
                                uniprotsnp.remove(tmp);
                            }
                        }
                    }
                    t=new StringTokenizer(work[6],";");
                    nt=t.countTokens();
                    for(int i=0;i<nt;i++){
                        String tmp=t.nextToken().toUpperCase();
                        if(genequery.containsKey("UNIPROT:"+tmp)){
                            if(allchr)chrsx.put(chr,"");
                            genes.put("UNIPROT:"+tmp,"UNIPROT:"+tmp);
                            genes.put(genename,"UNIPROT:"+tmp);
                            genequery.remove("UNIPROT:"+tmp);
                        }
                        if(hasuniprotsnp){
                            if(uniprotsnp.containsKey(tmp)){
                                if(allchr)chrsx.put(chr,"");
                                uniprotsnp.remove(tmp);
                            }
                        }
                    }
                }
                if(hasentrez){
                    t=new StringTokenizer(work[7],";");
                    nt=t.countTokens();
                    for(int i=0;i<nt;i++){
                        String tmp=t.nextToken().toUpperCase();
                        if(genequery.containsKey("ENTREZ:"+tmp)){
                            if(allchr)chrsx.put(chr,"");
                            genes.put(genename,"ENTREZ:"+tmp);
                            genequery.remove("ENTREZ:"+tmp);
                        }
                    }
                }
                StringBuilder sb2=new StringBuilder(work[3]);
                for(int i=4;i<work.length;i++)sb2.append("\t"+work[i]);
                String tmp=sb2.toString();
                dbNSFPgene.put(genename,compress(tmp));
                
            }
            in.close();
        }
        //gene not found in dbNSFP_gene
        if(genequery.size()>0) {
            genes.putAll(genequery);
            searchallchr=true;
        }
        if(uniprotsnp.size()>0)searchallchr=true;
        if(enspsnp.size()>0)searchallchr=true;
        if(hasrs)searchallchr=true;
        
        int nchr=chrs.length;
        if((!allchr)||(allchr&&!searchallchr)){
            nchr=chrsx.size();
            chrs=(String[])chrsx.keySet().toArray(new String[nchr]);
        }
        Arrays.sort(chrs);
        
        //begin search
        boolean printtitle=false;
        if(hg19)System.out.println("Searching based on hg19 coordinates:");
        else if(hg18)System.out.println("Searching based on hg18 coordinates:");
        else System.out.println("Searching based on hg38 coordinates:");
        for(int ii=0;ii<nchr;ii++){
            String chr=chrs[ii];
            System.out.println("\tSearching chr"+chr);
            infile2=classdir+"dbNSFP4.0b1a_variant.chr"+chr+".gz";
            f = new File(infile2);
            if (!f.isFile()){
                System.err.println("error: cannot find dbNSFP4.0b1a_variant.chr"+chr+".gz");
                continue;
            }
            in=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(infile2))));
            //int count=0;
            if(!printtitle) {
                if(isvcf&&preservevcf)out.print(title+"\t");
                if(allcol)out.println(in.readLine()+"\t"+genetitle);
                else{
                    String line=in.readLine()+"\t"+genetitle;
                    StringTokenizer t=new StringTokenizer(line,"\t");
                    int nt=t.countTokens();
                    String[] tmp=new String[nt];
                    for(int i=0;i<nt;i++)tmp[i]=t.nextToken();
                    out.print(tmp[cols[0]-1]);
                    for(int i=1;i<cols.length;i++){
                        out.print("\t"+tmp[cols[i]-1]);
                    }
                    out.println();
                }
                printtitle=true;
            }
            else in.readLine();
            String extra0="";
            String line;
            while((line=in.readLine())!=null){
                //count++;

                String[] temp=firstTokens(line, "\t", 18);
                String chr0=temp[0].toUpperCase();
                String pos=temp[1];
                String ref=temp[2].toUpperCase();
                String alt=temp[3].toUpperCase();
                String aaref=temp[4];
                String aaalt=temp[5];
                String rs=temp[6].toLowerCase();
                String hg19_chr=temp[7];
                String hg19_pos=temp[8];
                String hg18_chr=temp[9];
                String hg18pos=temp[10];
                String enspaapos=temp[11];
                String genename=temp[12].toUpperCase();
                String uniprotacc=temp[16].toUpperCase();
                String uniprotid=temp[17].toUpperCase();
                String uniprotaapos=enspaapos;
                String ensg=temp[13].toUpperCase();
                String enst=temp[14].toUpperCase();
                String ensp=temp[15].toUpperCase();
                String coor=chr0+"_"+pos;
                if(hg18) coor=hg18_chr+"_"+hg18pos;
                if(hg19) coor=hg19_chr+"_"+hg19_pos;
                
                boolean found=false;
                if(hasgene){
                    ArrayList names=new ArrayList();
                    StringTokenizer t3=new StringTokenizer(genename,"; ");
                    int nt3=t3.countTokens();
                    for(int i=0;i<nt3;i++) names.add(t3.nextToken());
                    if(hasuniprot) {
                        StringTokenizer t2=new StringTokenizer(uniprotid,";");
                        int nt2=t2.countTokens();
                        for(int i=0;i<nt2;i++){
                            uniprotid=t2.nextToken();
                            names.add("UNIPROT:"+uniprotid);
                        }
                        t2=new StringTokenizer(uniprotacc,";");
                        nt2=t2.countTokens();
                        for(int i=0;i<nt2;i++){
                            uniprotacc=t2.nextToken();
                            names.add("UNIPROT:"+uniprotacc);
                        }
                    }
                    if(hasens){
                        StringTokenizer t2=new StringTokenizer(ensg,";");
                        int nc1=t2.countTokens();
                        for(int i=0;i<nc1;i++) names.add("ENSEMBL:"+t2.nextToken());
                        t2=new StringTokenizer(enst,";");
                        nc1=t2.countTokens();
                        for(int i=0;i<nc1;i++) names.add("ENSEMBL:"+t2.nextToken());
                        t2=new StringTokenizer(ensp,";");
                        nc1=t2.countTokens();
                        for(int i=0;i<nc1;i++) names.add("ENSEMBL:"+t2.nextToken());
                    }
                    for(int i=0;i<names.size();i++){
                        String tmp=(String)names.get(i);
                        if(genes.containsKey(tmp)){
                            String tmp2=(String)genes.get(tmp);
                            foundgene.put(tmp2,null);
                            found=true;
                        }
                    }
                }
                String extra="";
                if(hasrs){
                    if(snprs.containsKey(rs)){
                        found=true;
                        snprs.put(rs,"");
                    }
                }
                if(hassnp){
                    String snp=coor+"_"+ref+"_"+alt;
                    String snp2=coor+"_"+ref+"_"+alt+"_"+aaref+"_"+aaalt;
                    
                    if(snps4.containsKey(snp))
                    {
                        found=true;
                        extra=(String)snps4.get(snp);
                        if(preservevcf)extra=decompress(extra);
                        extra=extra.trim();
                        if(preservevcf)snps4.put(snp, compress(extra));
                        else snps4.put(snp,"");
                    }
                    if(snps6.containsKey(snp2)){
                        found=true;
                        snps6.put(snp2, "");
                    }
                    if(poss.containsKey(coor)) {
                        found=true;
                        poss.put(coor," ");
                    }
                    if(hasuniprotsnp){
                        StringTokenizer t2=new StringTokenizer(uniprotid,";");
                        StringTokenizer t3=new StringTokenizer(uniprotacc,";");
                        StringTokenizer t4=new StringTokenizer(uniprotaapos,";");
                        int nt2=t2.countTokens();
                        int nt4=t4.countTokens();
                        String[] upos=new String[nt4];
                        for(int i=0;i<nt4;i++)upos[i]=t4.nextToken();
                        for(int i=0;i<nt2;i++){
                            uniprotid=t2.nextToken();
                            uniprotacc=t3.nextToken();
                            if(nt4==1)uniprotaapos=upos[0];
                            else uniprotaapos=upos[i];
                            String uniprot3="UNIPROT:"+uniprotid+":"+aaref+uniprotaapos+aaalt;
                            if(snpsp.containsKey(uniprot3)){
                                found=true;
                                snpsp.put(uniprot3, "");
                            }
                            String uniprot4="UNIPROT:"+uniprotacc+":"+aaref+uniprotaapos+aaalt;
                            if(snpsp.containsKey(uniprot4)){
                                found=true;
                                snpsp.put(uniprot4, "");
                            }
                        }
                    }
                    if(hasenspsnp){
                        StringTokenizer t2=new StringTokenizer(enst,";");
                        StringTokenizer t3=new StringTokenizer(ensp,";");
                        StringTokenizer t4=new StringTokenizer(enspaapos,";");
                        int nt2=t2.countTokens();
                        int nt4=t4.countTokens();
                        String[] epos=new String[nt4];
                        for(int i=0;i<nt4;i++)epos[i]=t4.nextToken();
                        for(int i=0;i<nt2;i++){
                            enst=t2.nextToken();
                            ensp=t3.nextToken();
                            if(nt4==1)enspaapos=epos[0];
                            else enspaapos=epos[i];
                            String ens3="ENSEMBL:"+enst+":"+aaref+enspaapos+aaalt;
                            if(snpsp.containsKey(ens3)){
                                found=true;
                                snpsp.put(ens3, "");
                            }
                            String ens4="ENSEMBL:"+ensp+":"+aaref+enspaapos+aaalt;
                            if(snpsp.containsKey(ens4)){
                                found=true;
                                snpsp.put(ens4, "");
                            }
                        }
                    }
                }
                
                if(found){
                    if(isvcf&&preservevcf){
                        if(extra.isEmpty())out.print(extra0+"\t");
                        else out.print(extra+"\t");
                    }
                    if(dbNSFPgene.containsKey(genename)){
                        String tmp=(String)dbNSFPgene.get(genename);
                        tmp=decompress(tmp);
                        if(allcol) out.println(line+"\t"+tmp);
                        else{
                            String linex=line+"\t"+tmp;
                            StringTokenizer tx=new StringTokenizer(linex,"\t");
                            int nt=tx.countTokens();
                            String[] tmpx=new String[nt];
                            for(int i=0;i<nt;i++)tmpx[i]=tx.nextToken();
                            out.print(tmpx[cols[0]-1]);
                            for(int i=1;i<cols.length;i++){
                                out.print("\t"+tmpx[cols[i]-1]);
                            }
                            out.println();
                        }
                    }
                    else {
                        String tmp=".";
                        for(int i=1;i<ngenecol-3;i++)tmp=tmp+"\t.";
                        if(allcol) out.println(line+"\t"+tmp);
                        else{
                            String linex=line+"\t"+tmp;
                            StringTokenizer tx=new StringTokenizer(linex,"\t");
                            int nt=tx.countTokens();
                            String[] tmpx=new String[nt];
                            for(int i=0;i<nt;i++)tmpx[i]=tx.nextToken();
                            out.print(tmpx[cols[0]-1]);
                            for(int i=1;i<cols.length;i++){
                                out.print("\t"+tmpx[cols[i]-1]);
                            }
                            out.println();
                        }
                    }
                }
                if(!extra.isEmpty())extra0=extra;
            }
            in.close();
        }
        out.close();

        //report unfound queries
        if(true){
            System.out.println();

            int nquery=snps4.size();
            String[] keys=(String[])snps4.keySet().toArray(new String[nquery]);
            int found=0;
            int unfound=0;
            for(int i=0;i<nquery;i++){
                String tmp=(String)snps4.get(keys[i]);
                if(preservevcf)tmp=decompress(tmp);
                if(tmp.startsWith(" ")) {
                    String tmp2=keys[i].replace("_", "\t");
                    out2.println(tmp2);
                    unfound++;
                }
                else found++;
            }
            nquery=snps6.size();
            keys=(String[])snps6.keySet().toArray(new String[nquery]);
            for(int i=0;i<nquery;i++){
                String tmp=(String)snps6.get(keys[i]);
                if(!tmp.isEmpty()) {
                    String tmp2=keys[i].replace("_", "\t");
                    out2.println(tmp2);
                    unfound++;
                }
                else found++;
            }
            nquery=snpsp.size();
            keys=(String[])snpsp.keySet().toArray(new String[nquery]);
            for(int i=0;i<nquery;i++){
                String tmp=(String)snpsp.get(keys[i]);
                if(!tmp.isEmpty()) {
                    String tmp2=keys[i];
                    out2.println(tmp2);
                    unfound++;
                }
                else found++;
            }
            nquery=snprs.size();
            keys=(String[])snprs.keySet().toArray(new String[nquery]);
            for(int i=0;i<nquery;i++){
                String tmp=(String)snprs.get(keys[i]);
                if(!tmp.isEmpty()) {
                    String tmp2=keys[i];
                    out2.println(tmp2);
                    unfound++;
                }
                else found++;
            }
            if(found>0)System.out.println(found+" SNP(s) are found. Written to "+outfile);
            if(unfound>0){
                System.out.println(unfound+" SNP(s) are not found. Written to "+outfile+".err");
                out2.println("Total "+unfound+" SNP(s) are not found.");
            }
            found=0;
            unfound=0;

            nquery=poss.size();
            keys=(String[])poss.keySet().toArray(new String[nquery]);
            for(int i=0;i<nquery;i++){
                String tmp=(String)poss.get(keys[i]);
                if(tmp.isEmpty()) {
                    String tmp2=keys[i].replace("_", "\t");
                    out2.println(tmp2);
                    unfound++;
                }
                else found++;
            }
            if(found>0)System.out.println(found+" position(s) are found. Written to "+outfile);
            if(unfound>0){
                System.out.println(unfound+" position(s) are not found. Written to "+outfile+".err");
                out2.println("Total "+unfound+" position(s) are not found.");
            }
            found=0;
            unfound=0;

            nquery=genes.size();
            String[] values=(String[])genes.values().toArray(new String[nquery]);
            for(int i=0;i<nquery;i++){
                if(!foundgene.containsKey(values[i]))unfoundgene.put(values[i],null);
            }
            found=foundgene.size();
            unfound=unfoundgene.size();
            keys=(String[])unfoundgene.keySet().toArray(new String[unfound]);
            for(int i=0;i<unfound;i++)out2.println(keys[i]);
            if(found>0)System.out.println(found+" gene(s) are found. Written to "+outfile);
            if(unfound>0){
                System.out.println(unfound+" gene(s) are not found. Written to "+outfile+".err");
                out2.println("Total "+unfound+" gene(s) are not found.");
            }
            out2.close();
        }
        
        //search splice region snp annotations
        
        if(searchsrsnp){
            if(!hg19&&!hg38){
                System.out.println("Error: to search dbscSNV, the coordinates must be based on hg19 or hg38.");
                System.exit(1);
            }
            out=new PrintWriter(new FileWriter(outfile2+".dbscSNV"));
            if(hg19)System.out.println("Searching scSNVs based on hg19");
            if(hg38)System.out.println("Searching scSNVs based on hg38");
            printtitle=false;
            for(int ii=0;ii<nchr;ii++){
                String chr=chrs[ii];
                System.out.println("\tSearching chr"+chr);
                infile2=classdir+"dbscSNV1.1.chr"+chr;
                f = new File(infile2);
                if (!f.isFile()){
                    System.err.println("error: dbscSNV1.1.chr"+chr+" not found.");
                    continue;
                }
                in=new BufferedReader(new FileReader(infile2));
                //int count=0;
                if(!printtitle) {
                    out.println(in.readLine()+"\t"+genetitle);
                    printtitle=true; 
                }
                else in.readLine();
                boolean ready=true;
                while(ready){
                    String line="";
                    try{
                            line = in.readLine();
                    }catch(java.io.IOException e){
                            //System.out.println(e);
                    }
                    if(line==null){

                        ready=false;
                        line="";
                        continue;
                    }
                    //count++;

                    StringTokenizer t=new StringTokenizer(line,"\t");
                    String chr0=t.nextToken().toUpperCase();
                    String pos=t.nextToken();
                    String ref=t.nextToken().toUpperCase();
                    String alt=t.nextToken().toUpperCase();
                    String chr38=t.nextToken().toUpperCase();
                    String pos38=t.nextToken();
                    for(int i=0;i<3;i++)t.nextToken();
                    String genename=t.nextToken();
                    int indx=genename.indexOf("(");
                    if(indx!=-1)genename=genename.substring(0,indx);
                    String coor=chr0+"_"+pos;
                    if(hg38)coor=chr38+"_"+pos38;
                    boolean found=false;
                    
                    if(hassnp){
                        String snp=coor+"_"+ref+"_"+alt;
                        
                        if(snps4.containsKey(snp))
                        {
                            found=true;
                        }
                        
                        if(poss.containsKey(coor)) {
                            found=true;
                        }
                        
                    }

                    if(found){
                        if(dbNSFPgene.containsKey(genename)){
                            String tmp=(String)dbNSFPgene.get(genename);
                            tmp=decompress(tmp);
                            out.println(line+"\t"+tmp);
                            
                        }
                        else {
                            String tmp=".";
                            for(int i=1;i<ngenecol-3;i++)tmp=tmp+"\t.";
                            out.println(line+"\t"+tmp);
                            
                        }
                    }
                }
                in.close();
            }
            out.close();
            String infile3="spidex_public_noncommercial_v1_0.tab.gz";
            String infile4="hg19_spidex.txt";
            if(hg19){
                f=new File(classdir+infile3);
                File f2=new File(classdir+infile4);
                if(f.isFile()){
                    out=new PrintWriter(new FileWriter(outfile2+".SPIDEX"));
                    System.out.println("Searching "+infile3);

                    in=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(classdir+infile3))));
                    out.println(in.readLine()+"\t"+genetitle);
                    boolean ready=true;
                    while(ready){
                        String line="";
                        try{
                                line = in.readLine();
                        }catch(java.io.IOException e){
                                //System.out.println(e);
                        }
                        if(line==null){

                            ready=false;
                            line="";
                            continue;
                        }
                        //count++;

                        StringTokenizer t=new StringTokenizer(line,"\t");
                        String chr=t.nextToken();
                        if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("chr")) chr=chr.substring(3);
                        if(chr.equalsIgnoreCase("MT"))chr="M";
                        if(chr.equalsIgnoreCase("XY"))chr="X";
                        String pos=t.nextToken();
                        String ref=t.nextToken();
                        String alt=t.nextToken();
                        for(int i=0;i<2;i++)t.nextToken();
                        String genename=t.nextToken();
                        String coor=chr+"_"+pos;
                        boolean found=false;

                        if(hassnp){
                            String snp=coor+"_"+ref+"_"+alt;

                            if(snps4.containsKey(snp))
                            {
                                found=true;
                            }

                            if(poss.containsKey(coor)) {
                                found=true;
                            }

                        }

                        if(found){
                            if(dbNSFPgene.containsKey(genename)){
                                String tmp=(String)dbNSFPgene.get(genename);
                                tmp=decompress(tmp);
                                out.println(line+"\t"+tmp);

                            }
                            else {
                                String tmp=".";
                                for(int i=1;i<ngenecol-3;i++)tmp=tmp+"\t.";
                                out.println(line+"\t"+tmp);

                            }
                        }
                    }
                    in.close();
                    out.close();
                    
                }
                else if(f2.isFile()){
                    out=new PrintWriter(new FileWriter(outfile2+".SPIDEX"));
                    System.out.println("Searching "+infile4);

                    in=new BufferedReader(new InputStreamReader((new FileInputStream(classdir+infile4))));
                    out.println(in.readLine());
                    boolean ready=true;
                    while(ready){
                        String line="";
                        try{
                                line = in.readLine();
                        }catch(java.io.IOException e){
                                //System.out.println(e);
                        }
                        if(line==null){

                            ready=false;
                            line="";
                            continue;
                        }
                        //count++;

                        StringTokenizer t=new StringTokenizer(line,"\t");
                        String chr=t.nextToken();
                        if(chr.length()>3&&chr.substring(0, 3).equalsIgnoreCase("chr")) chr=chr.substring(3);
                        if(chr.equalsIgnoreCase("MT"))chr="M";
                        if(chr.equalsIgnoreCase("XY"))chr="X";
                        String pos=t.nextToken();
                        String end=t.nextToken();
                        String ref=t.nextToken();
                        String alt=t.nextToken();
                        String coor=chr+"_"+pos;
                        boolean found=false;

                        if(hassnp){
                            String snp=coor+"_"+ref+"_"+alt;

                            if(snps4.containsKey(snp))
                            {
                                found=true;
                            }

                            if(poss.containsKey(coor)) {
                                found=true;
                            }

                        }

                        if(found){
                            out.println(line);
                        }
                    }
                    in.close();
                    out.close();
                }
                else{
                    System.out.println("spidex_public_noncommercial_v1_0.tab.gz or hg19_spidex.txt is not found.");
                }
                    
            }
            
        }
    }

    private  static String compress(String str) throws IOException {
        if (str == null || str.length() == 0) {
            return str;
        }
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        GZIPOutputStream gzip = new GZIPOutputStream(out);
        gzip.write(str.getBytes());
        gzip.close();
        String outStr = out.toString("ISO-8859-1");
        return outStr;
     }
    
    private  static String decompress(String str) throws IOException {
        if (str == null || str.length() == 0) {
            return str;
        }
        GZIPInputStream gis = new GZIPInputStream(new ByteArrayInputStream(str.getBytes("ISO-8859-1")));
        BufferedReader bf = new BufferedReader(new InputStreamReader(gis, "ISO-8859-1"));
        String outStr = "";
        String line;
        while ((line=bf.readLine())!=null) {
          outStr += line;
        }
        return outStr;
     }
    private static String[] firstTokens(String s, String delimiter, int num){  
    if (s == null) {  
        return null;  
    }  
    int stringLength = s.length();
    int delimiterLength;  
    if (delimiter == null || (delimiterLength = delimiter.length()) == 0){  
        return new String[] {s};  
    }  

    int start;  
    int end;  
  
  
    String[] result = new String[num];  
    Arrays.fill(result, "");
    
    int count = 0;  
    start = 0;  
    while((end = s.indexOf(delimiter, start)) != -1&&count<num){  
        result[count] = (s.substring(start, end));  
        count++;  
        start = end + delimiterLength;  
    }  
    if(count<num){
        end = stringLength;  
        result[count] = s.substring(start, end);
    }
          
  
    return (result);  
}
    private static boolean isInt(String s) {
try {
Integer.parseInt(s);
}
catch (NumberFormatException nfe) {
return false;
}
return true;
}
}
