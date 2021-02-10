REG := 100225593120.dkr.ecr.us-east-1.amazonaws.com
DOCKER_IMAGE_TAG := latest

registry-docker-login:
ifneq ($(shell echo ${REG} | egrep "ecr\..+\.amazonaws\.com"),)
	@$(eval DOCKER_LOGIN_CMD=aws)
ifneq (${AWS_PROFILE},)
	@$(eval DOCKER_LOGIN_CMD=${DOCKER_LOGIN_CMD} --profile ${AWS_PROFILE})
endif
	@$(eval DOCKER_LOGIN_CMD=${DOCKER_LOGIN_CMD} ecr get-login-password | docker login -u AWS --password-stdin https://${REG})
	${DOCKER_LOGIN_CMD}
endif

build: registry-docker-login
	docker build -t ${REG}/agr_preprocess_run:${DOCKER_IMAGE_TAG} --build-arg REG=${REG} .

run: build
	REG=${REG} DOCKER_IMAGE_TAG=${DOCKER_IMAGE_TAG} docker-compose up agr_preprocess

bash:
	REG=${REG} DOCKER_IMAGE_TAG=${DOCKER_IMAGE_TAG} docker-compose up agr_preprocess bash
