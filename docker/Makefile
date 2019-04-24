TAG=latest  # tag to give docker container on build
VERSION=master  # version of hscPipe to build

all:	hscPipe

hscPipe:	Dockerfile
	docker build -t paprice/hscpipe:$(TAG) -f Dockerfile --build-arg HSCPIPE_VERSION=$(VERSION) .

# For forcing a completely new build
hscPipe-force:	Dockerfile
	docker build -t paprice/hscpipe:$(TAG) -f Dockerfile --build-arg HSCPIPE_VERSION=$(VERSION) --build-arg FORCE=$$(date +%s) .

.PHONY: hscPipe hscPipe-force
