brouwer.txt: brouwer.tmp
	sage tmp_to_db.py
	rm -f brouwer.tmp

brouwer.tmp:
	./website_to_tmp.sh

clean:
	rm -f brouwer.txt
	rm -f brouwer.sobj
	rm -f brouwer.tmp

