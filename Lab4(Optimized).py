import pprint

condon_tablet={
'UUU':'F','CUU':'L','AUU':'I','GUU':'V',
'UUC':'F','CUC':'L','AUC':'I','GUC':'V',
'UUA':'L','CUA':'L','AUA':'I','GUA':'V',
'UUG':'L','CUG':'L','AUG':'M','GUG':'V',
'UCU':'S','CCU':'P','ACU':'T','GCU':'A',
'UCC':'S','CCC':'P','ACC':'T','GCC':'A',
'UCA':'S','CCA':'P','ACA':'T','GCA':'A',
'UCG':'S','CCG':'P','ACG':'T','GCG':'A',
'UAU':'Y','CAU':'H','AAU':'N','GAU':'D',
'UAC':'Y','CAC':'H','AAC':'N','GAC':'D',
'UAA':'|','CAA':'Q','AAA':'K','GAA':'E',
'UAG':'|','CAG':'Q','AAG':'K','GAG':'E',
'UGU':'C','CGU':'R','AGU':'S','GGU':'G',
'UGC':'C','CGC':'R','AGC':'S','GGC':'G',
'UGA':'|','CGA':'R','AGA':'R','GGA':'G',
'UGG':'W','CGG':'R','AGG':'R','GGG':'G'
}

com_base={
    'A':'U',
    'T':'A',
    'C':'G',
    "G":'C'
}

def read_str(str):
    condons=[]
    for i in range(0,len(str),3):
        condons.append((str[i:i+3].upper()))
    return condons

def transcript(DNA=str):
    mRNA=''
    for ji in DNA.upper():
        mRNA+=com_base[ji]
    return mRNA

def translate(str):
    condons=read_str(str)
    AAs=''
    for condon in condons:
        try:
            AA=condon_tablet[condon]
        except:
            AA='*'
        AAs+=AA
    return AAs




def file_process(abs_path,outSOURCE='1'):
    data={
        'REFERENCE':{},
        'LOCUS':'',
        'DEFINITION':'',
        'ACCESSION':'',
        'VERSION':'',
        'KEYWORDS':'',
        'COMMENT':'',
        'ORIGIN':'',
        'SOURCE':'',
        'PRIMARY' :'', 
        'FEATURES' :'',
        'ORGANISM':''
    }
    file=open(abs_path)
    lines=file.readlines()
    #REFS need sub_FLAG
    main_index=['COMMENT','PRIMARY','FEATURES',"ORGANISM"]
    single_index=['LOCUS' , 'DEFINITION' , 'ACCESSION' , 'VERSION' , 'KEYWORDS','SOURCE']
    upper_index=main_index+single_index+["ORIGIN"]+["REFERENCE"]
    detail=''
    #No REF and Origin
    i=-1
    main_FLAG=''
    sub_FLAG=''
    ref_index=0
    ref_components=['AUTHORS','TITLE','JOURNAL','PUBMED','REMARK']
    for line in lines:
        ls=line.strip().split(' ')
        #可选择是否保留原有文件格式,此处不保留
        head=ls[0]
        details=''
        if head in upper_index+ref_components:
            #LOCUS and so on.
            if head in ['LOCUS' , 'DEFINITION' , 'ACCESSION' , 'VERSION' , 'KEYWORDS']:
                for detail in ls[1:]:
                    details+=detail+' '
                data[head]+=details

            elif head in ['COMMENT' , 'SOURCE' , 'PRIMARY' , 'FEATURES' , 'ORIGIN' , 'ORGANISM']:
                main_FLAG=head
                for detail in ls[1:]:
                    details+=detail+' '
                data[head]+=details
            #To create or update a REFERENCE
            elif head=='REFERENCE':
                main_FLAG=head
                ref_index+=1
                exec('Reference'+str(ref_index)+"={}")
                for i in ref_components:
                    eval('Reference'+str(ref_index))[i]=''
            elif head in ref_components:
                sub_FLAG=head
                for detail in ls[1:]:
                    details+=detail+' '
                eval('Reference'+str(ref_index))[head]+=details         
        else:
            if main_FLAG=='ORIGIN':
                for detail in ls[1:]:
                    details+=detail+''
                data[main_FLAG]+=details
            elif main_FLAG=='REFERENCE' :
                for detail in ls:
                    details+=detail+' '
                eval('Reference'+str(ref_index))[sub_FLAG]+=details
                #！
            else:
                data[main_FLAG]+=' '+line.strip()+' '
    for i in range(1,ref_index+1):
        data['REFERENCE'][i]=(eval('Reference'+str(i)))
    pprint.pprint(data)



    if outSOURCE=='1':
        file=open('Source of '+data["ACCESSION"].strip()+'.FASTA','w')
        file.write(f">{data['ACCESSION']}\n")
        for j in range(0,len(data['ORIGIN']),60):
            file.write(data['ORIGIN'][j:j+60]+'\n')
        file.close()

    return data
