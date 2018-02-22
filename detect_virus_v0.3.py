import os
import sys
import time
import threading
import thread
from math import log
from time import sleep,ctime
import multiprocessing
from optparse import OptionParser
from detect_virus import *
def multisample(tag,bam_dic,virus_list,config,file_list=[]):
    l=virus_list
    outfile=open(tag+'.txt','w')
    if config=='N':
        filelist=os.listdir(bam_dic)
        for filename in filelist:
            if filename.startswith(tag) and filename.endswith('bam'):
                for line in l:
                    cmd='samtools view '+bam_dic+'/'+filename+' '+line+' > test.txt'
                    os.system(cmd)
                    for tmp in open('test.txt'):
                        if not tmp.startswith('@'):
                            newline=tmp.rstrip().split('\t')
                            newline.append(filename)
                            if (newline[6] !='=' or 'S' in newline[5]) and newline[6] not in l:
                                outfile.write('\t'.join(newline)+'\n')
    else:
        
        for filename in file_list:
            for line in l:
                cmd='samtools view '+bam_dic+'/'+filename+' '+line+' > test.txt'
                os.system(cmd)
                for tmp in open('test.txt'):
                    if not tmp.startswith('@'):
                        newline=tmp.rstrip().split('\t')
                        newline.append(filename)
                        if (newline[6] !='=' or 'S' in newline[5]) and newline[6] not in l:
                            outfile.write('\t'.join(newline)+'\n')
    outfile.close()
def compare(l,dic,key):
    tag=0
    for i in range(len(dic[key])):
        if l[0]<=dic[key][i][0]<=l[1]<=dic[key][i][1]:
            dic[key][i][0]=l[0]
            dic[key][i][-1]=dic[key][i][-1]+1
            dic[key][i][-2].append(l[-2][0])
            dic[key][i][2]=min(dic[key][i][2],l[2])
            dic[key][i][3]=max(dic[key][i][3],l[3])
            tag=1
            break
        elif dic[key][i][0]<=l[0]<=dic[key][i][1]<=l[1]:
            dic[key][i][1]=l[1]
            dic[key][i][-1]=dic[key][i][-1]+1
            dic[key][i][-2].append(l[-2][0])
            dic[key][i][2]=min(dic[key][i][2],l[2])
            dic[key][i][3]=max(dic[key][i][3],l[3])
            tag=1
            break
        elif l[0]<=dic[key][i][0]<=dic[key][i][1]<=l[1]:
            dic[key][i][0]=l[0]
            dic[key][i][1]=l[1]
            dic[key][i][-1]=dic[key][i][-1]+1
            dic[key][i][-2].append(l[-2][0])
            dic[key][i][2]=min(dic[key][i][2],l[2])
            dic[key][i][3]=max(dic[key][i][3],l[3])
            tag=1
            break
        elif dic[key][i][0]<=l[0]<=l[1]<=dic[key][i][1]:
            dic[key][i][-1]=dic[key][i][-1]+1
            dic[key][i][-2].append(l[-2][0])
            dic[key][i][2]=min(dic[key][i][2],l[2])
            dic[key][i][3]=max(dic[key][i][3],l[3])
            tag=1
            break
        else:
            continue
    if tag==1:
        return [dic,True]
    else:
        return [dic,False]
    
def read_file(filename,n=300):
    dic={}
    for line in open(filename):
        if not line.startswith('@'):
            newline=line.rstrip().split('\t')
            if newline[6] !='=':
                if newline[6] not in dic:
                    dic[newline[6]]=[[int(newline[7])-n,int(newline[7])+n,max(int(newline[3])-n,1),int(newline[3])+n,newline[0],newline[2],[newline[-1]],1]]                 
                else:
                    result=compare([int(newline[7])-n,int(newline[7])+n,max(int(newline[3])-n,1),int(newline[3])+n,newline[0],newline[2],[newline[-1]],1],dic,newline[6])
                    if result[-1] is True:
                        dic=result[0]
                    else:
                        dic[newline[6]].append([int(newline[7])-n,int(newline[7])+n,max(int(newline[3])-n,1),int(newline[3])+n,newline[0],newline[2],[newline[-1]],1])
            else:
                continue
        else:
            continue
    return dic

def read_soft_clip(filename,split_length):
    soft_list=[]
    for line in open(filename):
        newline=line.rstrip().split('\t')
        S_num=newline[5].count('S')
        if S_num==1:
            S_tag=newline[5].split('S')
            if S_tag[-1]=='':
                S_len=S_tag[0].split('M')[-1]
                if int(S_len)>=split_length:
                    soft_list.append(line)
                else:
                    continue
            else:
                if int(S_tag[0])>=split_length:
                    soft_list.append(line)
                else:
                    continue
        elif S_num==2:
            S_tag=newline[5].split('S')
            S_len=S_tag[1].split('M')[-1]
            if int(S_len)>=split_length or int(S_tag[0])>=split_length:
                soft_list.append(line)
            else:
                continue
        else:
            print newline[5]
    
    return soft_list

def extract_clip(dic,dic3,filename,bam_dic,tag,mini_support,split_length,n=300):
    result=[]
    err_log=open('error.txt','w')
    outwrong=open('wrong_SV.txt','w')
    #filelist=os.listdir(bam_dic)
    for key in dic:
        for line in dic[key]:
            if line[-1]>mini_support:
                
                line[-2]=list(set(line[-2]))
                tmp_line=[[key,str((line[0]+line[1])/2),line[5],str((line[2]+line[3])/2)]]
                #print [key,line,tmp_line]
                result_score=sys_exe(tmp_line,dic3,bam_dic,line[-2],n,split_length,0)
                print result_score
                if result_score[0]<result_score[1]:
                    outfile=open('soft_clipping_'+key+'_'+str(line[0])+'_'+str(line[1])+'.txt','w')
                    #for t in filelist:
                    #    if t.startswith(tag) and t.endswith('bam') and t not in line[-2]:
                    #        line[-2].append(t)
                    for tmp_file in line[-2]:                     
                        cmd='samtools view '+bam_dic+'/'+tmp_file+' '+key+':'+str(line[0])+'-'+str(line[1])+' > tmt.txt'
                        #print cmd
                        os.system(cmd)
                        for tmp_line in open('tmt.txt'):
                            if not tmp_line.startswith('@'):
                                new_tmp_line=tmp_line.rstrip().split('\t')
                                if 'S' in new_tmp_line[5] and line[0]<=int(new_tmp_line[3])<=line[1]:
                                    outfile.write(tmp_line.rstrip()+'\t'+tmp_file+'\n')
                                else:
                                    continue
                            else:
                                continue
                        for tmp_line in open(filename):
                            if not tmp_line.startswith('@'):
                                new_tmp_line=tmp_line.rstrip().split('\t')
                                if 'S' in new_tmp_line[5] and line[2]<=int(new_tmp_line[3])<=line[3] and new_tmp_line[2]==line[-3]:
                                    outfile.write(tmp_line.rstrip()+'\t'+tmp_file+'\n')
                                else:
                                    continue
                            else:
                                continue
                    outfile.close()
                    soft_list=read_soft_clip('soft_clipping_'+key+'_'+str(line[0])+'_'+str(line[1])+'.txt',split_length)
                    #print [key,soft_list]
                    result.append([key,line[0],line[1],line[5],line[2],line[3],soft_list])
                elif result_score==[0,0]:
                    outwrong.write('\t'.join(tmp_line[0])+'\terror\n')
                else:
                    outwrong.write('\t'.join(tmp_line[0])+'\n')
            else:
                continue
    err_log.close()
    outwrong.close()
    return result
def generate_fasta(result,dic,nloop,split_length):
    outfile_new=open('virus_integration_'+str(nloop+1)+'.txt','w')
    #outfile_new.write('chrome\tpois\tvirus\tpois\n')
    for line in result:
        #print line
        outfile=open('ref_'+line[0]+'_'+str(nloop+1)+'.fa','w')
        outfile1=open('seq_'+line[0]+'_'+str(nloop+1)+'.fa','w')
        str_human=dic[line[0]][line[1]-1:line[2]].upper()
        str_virus=dic[line[3]][line[4]-1:line[5]].upper()
        ref=str_human+str_virus
        outfile.write('>ref\n')
        outfile.write(ref)
        for new_tmp in line[-1]:
            tmp=new_tmp.rstrip().split('\t')
            outfile1.write('>'+tmp[0]+'_'+tmp[-1]+'\n')
            outfile1.write(tmp[9].upper()+'\n')
        outfile.close()
        outfile1.close()
        cmd='pairwise_align -inv '+'ref_'+line[0]+'_'+str(nloop+1)+'.fa'+' '+'seq_'+line[0]+'_'+str(nloop+1)+'.fa'+' >align_'+line[0]+'_'+str(nloop+1)+'.txt'
        #print cmd
        os.system(cmd)
        read_result('align_'+line[0]+'_'+str(nloop+1)+'.txt',line,outfile_new,split_length)
    outfile_new.close()
def reference(ref_name):
    ref_dic={}
    chr_name=''
    for line in open(ref_name):
        newline=line.rstrip()
        if newline.startswith('>'):
            if chr_name!='':
                ref_dic[chr_name]=tmp_str
            chr_name=newline.split('>')[1]
            tmp_str=''
        else:
            tmp_str=tmp_str+newline
    ref_dic[chr_name]=tmp_str
    return ref_dic
def read_result(filename,pois,outfile,split_length):
    result=[]
    for line in open(filename):
        if 'first' in line and 'EXCISED REGION' in line:
            newline=line.rstrip()
            tmp1_left=int(newline.split('=>  ')[1].split(' EXCISED REGION ')[0].split('[')[1].split(']')[0].split(',')[0])
            tmp1_right=int(newline.split('=>  ')[1].split(' EXCISED REGION ')[0].split('[')[1].split(']')[0].split(',')[1])
            tmp2_left=int(newline.split('=>  ')[1].split(' EXCISED REGION ')[1].split('[')[1].split(']')[0].split(',')[0])
            tmp2_right=int(newline.split('=>  ')[1].split(' EXCISED REGION ')[1].split('[')[1].split(']')[0].split(',')[1])
            if int(score)>=37 and abs(int(tmp1_right)-int(tmp1_left))>=split_length and abs(int(tmp2_right)-int(tmp2_left))>=split_length:
                if max(tmp1_right,tmp1_left)<pois[2]-pois[1]:
                    if min(tmp2_right,tmp2_left)>=pois[2]-pois[1]:
                        result.append([tmp1_left,tmp1_right,tmp2_left,tmp2_right])
                else:
                    if max(tmp2_right,tmp2_left)<pois[2]-pois[1]:
                        result.append([tmp1_left,tmp1_right,tmp2_left,tmp2_right])
        elif line.startswith('Aligned:'):
            score=line.split('Aligned:        ')[1].split('    ')[0]
        else:
            continue
    left=[]
    right=[]
    for line in result:
        temp1=max(int(line[0]),int(line[1]))
        temp2=max(int(line[2]),int(line[3]))
        if temp1>pois[2]-pois[1]:
            right.append(temp1)
            left.append(temp2)
        else:
            left.append(temp1)
            right.append(temp2)
    dic_left={}
    dic_right={}
    for tmp_left in sorted(left,reverse=True):
        if tmp_left in dic_left:
            dic_left[tmp_left]=dic_left[tmp_left]+1
        else:
            key=dic_left.keys()
            if key==[]:
                dic_left[tmp_left]=1
            else:
                for tmp_key in key:
                    if abs(tmp_left-tmp_key)<=20:
                        dic_left[tmp_key]=dic_left[tmp_key]+1
                    else:
                        dic_left[tmp_left]=1
    for tmp_right in sorted(right,reverse=True):
        if tmp_right in dic_right:
            dic_right[tmp_right]=dic_right[tmp_right]+1
        else:
            key=dic_right.keys()
            if key==[]:
                dic_right[tmp_right]=1
            else:
                tmp_tag=0
                for tmp_key in key:
                    if abs(tmp_right-tmp_key)<=20:
                        dic_right[tmp_key]=dic_right[tmp_key]+1
                        tmp_tag=1
                        break
                    else:
                        continue
                if tmp_tag==0:
                    dic_right[tmp_right]=1
                        
    sort_left=sorted(dic_left.iteritems(),key=lambda x:(x[1],x[0]),reverse=True)
    sort_right=sorted(dic_right.iteritems(),key=lambda x:(x[1],x[0]),reverse=True)
    print [sort_left,pois[:6]]
    print [sort_right,pois[:6]]
    if sort_left==[] or sort_right==[]:
        print 'no integration'
    else:
        if sort_left[0][1]>=2 and sort_right[0][1]>=2:
            integration_pois=sort_left[0][0]+pois[1]-1
            virus_pois=sort_right[0][0]-(pois[2]-pois[1])+pois[4]-1
        
            print[integration_pois,virus_pois]
            outfile.write(pois[0]+'\t'+str(integration_pois)+'\t'+pois[3]+'\t'+str(max(virus_pois,1))+'\n')
        else:
            print 'not integration'
def splice(L,n):
    tmp_L=[[] for i in range(n)]
    count=0
    for line in L:
        m=count%n
        tmp_L[m].append(line)
        count=count+1
    return tmp_L
def read_config(filename):
    result=[]
    for line in open(filename):
        result.append(line.rstrip())
    return result
def main():
    usage = """%prog -h <multi tag of bam files> -d <mulit bamfile dictionary> -r <human reference> -v <virus reference> -t <multi thread> -o <output>
detect_viurs_v0.1.0
Author: Yuchao Xia	
Description: detect integration of virus. 
	"""
    print 'starting at:',ctime()
    time_total=time.time()
    parser = OptionParser(usage)
	
    parser.add_option("-i", "--inFile", dest="head", help="head of multi bam files",metavar="STR")
    parser.add_option("-d", "--dic", dest="bam_dic", help="dictionary of mutii bam files",metavar="FILE")
    parser.add_option("-t", "--thread", dest="thread", default='1',help="mulitithread",metavar='INT')
    parser.add_option("-r", "--ref", dest="ref",help="reference of human",metavar="FILE")
    parser.add_option("-v", "--virus", dest="virus",help="reference of virus")
    parser.add_option("-o", "--output", dest="output", default='integration.txt',help="outfile of integration file")
    parser.add_option("-p", "--pair", dest="pair", default='2',help="the minimum number of supported reads")
    parser.add_option("-s", "--split", dest="split", default='15',help="the minimum length of split reads")
    parser.add_option("-w", "--window", dest="window", default='300',help="the length of window")
    parser.add_option("-c", "--config", dest="config",default='',help="the name of bam file")
    (ops, args) = parser.parse_args()
    if ops.head is None or ops.bam_dic is None or ops.ref is None:
        parser.print_help()
    else:
        print ("minisupport reads: ", ops.pair)
        print ('minisplit reads:',ops.split)
        print('window:',ops.window)
        tag=ops.head
        dic3=reference(ops.ref)
        if ops.virus is None:
            virus_list=[]
            key=dic3.keys()
            for tmp in key:
                if not tmp.startswith('chr'):
                    virus_list.append(tmp)
                else:
                    continue
        else:
            dic4=reference(ops.virus)
            virus_list=dic4.keys()
        if ops.config!='':
            config_file=read_config(ops.config)
            multisample(tag,ops.bam_dic,virus_list,'Y',config_file)
        else:
            multisample(tag,ops.bam_dic,virus_list,'N')
        dic=read_file(tag+'.txt')
        result=extract_clip(dic,dic3,tag+'.txt',ops.bam_dic,tag,int(ops.pair),int(ops.split),int(ops.window))
        result=result[:30]
        sv_ctx_splice=splice(result,int(ops.thread))
        for line in sv_ctx_splice:
            print(len(line))
            
        process=[]
        nloops=range(len(sv_ctx_splice))
        for i in nloops:
            t=multiprocessing.Process(target= generate_fasta,args=(sv_ctx_splice[i],dic3,i,int(ops.split)))
            process.append(t)
        for i in nloops:
            process[i].start()
        for i in nloops:
            process[i].join()
        outfile=open(ops.output,'w')
        outfile.write('#chrome\tpois\tvirus\tpois\n')
        outfile.close()
        cmd='cat virus_integration_*.txt >> '+ops.output
        os.system(cmd)
        cmd1='rm soft_clipping*.txt'
        cmd2='rm ref*.fa'
        cmd3='rm seq*.fa'
        cmd4='rm align*.txt'
        cmd5='rm virus_integration_*.txt'
        cmd6='tmt.txt'
        cmd7='test.txt'
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
        #os.system(cmd6)
        #os.system(cmd7)
        sort_l=[]
        for line in open(ops.output):
            if not line.startswith('#'):
                newline=line.rstrip().split('\t')
                sort_l.append(newline)
            else:
                continue
        sort_l_new=sorted(sort_l,key=lambda d:(d[0],int(d[1])))
        outfile=open(ops.output,'w')
        outfile.write('#chrome\tpois\tvirus\tpois\n')
        for line in sort_l_new:
            outfile.write('\t'.join(line)+'\n')
        outfile.close()
    time_end_total=time.time()
    elapsed_total=time_end_total-time_total
    print ("Time taken: ", elapsed_total, "seconds.")         
    print 'all done at:',ctime()
    #generate_fasta(result,dic3,nloop)
main()
                
#python soft_clip.py 213.txt /data1/XiLab/qfxing/upload_HCC/output1                        
                
            
