build:
	docker build -t agrdocker/agr_preprocess_run:latest .

run: build
	docker-compose up agr_preprocess

bash:
	docker-compose up agr_preprocess bash
