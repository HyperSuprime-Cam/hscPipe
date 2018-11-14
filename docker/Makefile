TAG=latest
VERSION=6.7

hscPipe:	Dockerfile
	docker build -t paprice/hscpipe:$(TAG) -f Dockerfile --build-arg HSCPIPE_VERSION=$(VERSION) --build-arg FORCE=$$(date +%s) .


all:	hscPipe
.PHONY: hscPipe
