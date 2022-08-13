# !/bin/bash
# 用法：bash wgetSRA.sh list
# 其中list为存放SRR/ERR/DRR编号的文件，一行一个ID
# 依赖：SRA Toolkit 和 mwget，下载，配置，添加至环境变量
# 依赖检测：srapath -h 和 vdb-validate -h 和 mwget -h


if [ -a $1 ];then

for i in $(cat $1);do

if [ $(grep -c $i success.fin.txt) -eq "1" ]; then
echo -e "\033[32m${i}已经下载成功，详情请查看success.fin.txt\033[0m"
sed -i "/${i}/d" $1
else
mwget -n 8 -f ${i} $(srapath ${i})
vdb-validate ${i} > ${i}.log 2>&1
if [ $(tail -1 ${i}.log | awk -vFS=" " '{print $NF}') == "consistent" ];then
echo -e "\033[32m${i}校验成功，详情请查看success.fin.txt\033[0m"
sed -i "/${i}/d" $1
echo -e "$(date)\t$i" >> success.temp.txt
awk -vFS="\t" '!a[$2]++' success.temp.txt > success.fin.txt
else
echo -e "\033[33m${i}校验失败，重新拉起下载\033[0m"
rm ${i}
rm ${i}.mg!
mwget -n 8 -f ${i} $(srapath ${i})
vdb-validate ${i} >${i}.log 2>&1
if [ $(tail -1 ${i}.log | awk -vFS=" " '{print $NF}') == "consistent" ];then
echo -e "\033[32m${i}重新下载后校验成功，详情请查看success.fin.txt\033[0m"
sed -i "/${i}/d" $1
echo -e "$(date)\t$i" >> success.temp.txt
awk -vFS="\t" '!a[$2]++' success.temp.txt > success.fin.txt
else
echo -e "\033[31m${i}再次校验失败，取消下载，详情请查看fail.fin.txt\033[0m"
echo -e "$(date)\t$i" >> fail.temp.txt
awk -vFS="\t" '!a[$2]++' fail.temp.txt > fail.fin.txt
fi
fi
rm ${i}.log
fi
done

else
echo "error：找不到名为$1的文件，请修改存放SRR_ID的文件名后重试"
exit 1

fi
