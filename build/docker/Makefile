CUR_DATE=`date +'%m_%d_%Y'`

all: build_docker_image login_to_dockerhub push_to_dockerhub

build_docker_image:
	docker build --no-cache -t aka1spanky/base_image:$(CUR_DATE) .

login_to_dockerhub:
	docker login --username=aka1spanky

push_to_dockerhub:
	docker push aka1spanky/base_image
