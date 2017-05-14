Commands to run:

./rank_percent.py star logs/WABI_2017/star/bounded/ > logs/star_bounded
./newf1.py -k 10 -i logs/WABI_2017/star/bounded/ -n star >> logs/star_bounded
cat logs/star_bounded | python latex_tables.py | sublime
