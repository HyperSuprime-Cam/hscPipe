TAG=latest
VERSION=7.0

all:	hscPipe

hscPipe:	Dockerfile
	docker build -t paprice/hscpipe:$(TAG) -f Dockerfile --build-arg HSCPIPE_VERSION=$(VERSION) .

# For forcing a completely new build
hscPipe-force:	Dockerfile
	docker build -t paprice/hscpipe:$(TAG) -f Dockerfile --build-arg HSCPIPE_VERSION=$(VERSION) --build-arg FORCE=$$(date +%s) .

.PHONY: hscPipe hscPipe-force
