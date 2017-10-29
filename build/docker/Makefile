CUR_DATE=`date +'base_hub_%m_%d_%Y'`

docker_build_jupyterhub_image:
	docker build --no-cache -t aka1spanky/fractal_jupyterhub:$(CUR_DATE) .

login_to_dockerhub:
	docker login --username=aka1spanky

push_jupyterhub_image_to_dockerhub:
	docker push aka1spanky/fractal_jupyterhub:$(CUR_DATE)

