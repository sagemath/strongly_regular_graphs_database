database:
	./website_to_tmp.sh
	python tmp_to_db.py
	rm -f brouwer.tmp
clean:
	rm -f brouwer.txt
	rm -f brouwer.json
	rm -f brouwer.tmp

