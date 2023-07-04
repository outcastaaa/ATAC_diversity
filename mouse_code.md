```bash
cd ~/xuruizhi/brain/brain/trim/mouse/ 
mkdir -p ~/xuruizhi/brain/brain/alignment_new/mouse
# 写双端批量
cat 111.sh
#!/usr/bin bash
bowtie2  -p 48 -x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -1 {}_1_val_1.fq.gz -2 {}_2_val_2.fq.gz \
2> ~/xuruizhi/brain/brain/alignment_new/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment_new/mouse/{}.sam

cat pair.list  | while read id; do sed "s/{}/${id}/g" 111.sh > ${id}.sh; done
# Job <8601500> is submitted to queue <mpi>.
# Job <8601501> is submitted to queue <mpi>.
# Job <8601502> is submitted to queue <mpi>.
# Job <8601503> is submitted to queue <mpi>.
# Job <8601504> is submitted to queue <mpi>.
# Job <8601505> is submitted to queue <mpi>.
# Job <8601506> is submitted to queue <mpi>.
# Job <8601507> is submitted to queue <mpi>.
# Job <8601508> is submitted to queue <mpi>.
# Job <8601509> is submitted to queue <mpi>.
# Job <8601510> is submitted to queue <mpi>.
# Job <8601511> is submitted to queue <mpi>.
# Job <8601512> is submitted to queue <mpi>.
# Job <8601513> is submitted to queue <mpi>.
# Job <8601514> is submitted to queue <mpi>.
# Job <8601515> is submitted to queue <mpi>.
# Job <8601516> is submitted to queue <mpi>.


# 写单端批量
cat 222.sh
#!/usr/bin bash
bowtie2  -p 48 -x  ~/xuruizhi/brain/brain/genome/mouse/mm10 \
--very-sensitive -X 2000 -U {}_trimmed.fq.gz \
2> ~/xuruizhi/brain/brain/alignment_new/mouse/{}.summary \
-S ~/xuruizhi/brain/brain/alignment_new/mouse/{}.sam

cat single.list  | while read id; do sed "s/{}/${id}/g" 222.sh > ${id}.sh; done

cat single.list | while read id
do
  bsub -q mpi -n 48 -o ~/xuruizhi/brain/brain/alignment_new/mouse/ "bash ${id}.sh"
done
# Job <8602391> is submitted to queue <mpi>.
# Job <8602392> is submitted to queue <mpi>.
# Job <8602393> is submitted to queue <mpi>.
```